// EM algorithm for MNHMMs

#include "mnhmm.h"
#include "logsumexp.h"
#include "sum_to_zero.h"
#include <nloptrAPI.h>
#include <chrono>

double mnhmm::objective_pi(const arma::vec& x, arma::vec& grad) {
  
  mstep_iter++;
  double value = 0;
  eta_pi(current_d) = arma::mat(x.memptr(), S - 1, X_pi.n_rows);
  gamma_pi(current_d) = sum_to_zero(eta_pi(current_d), Qs);
  grad.zeros();
  arma::mat tQs = Qs.t();
  arma::uvec idx(S);
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i, current_d);
    }
    const arma::vec& counts = E_pi(current_d).col(i);
    idx = arma::find(counts);
    double val = arma::dot(counts.rows(idx), log_pi(current_d).rows(idx));
    if (!std::isfinite(val)) {
      grad.zeros();
      return maxval;
    }
    value -= val;
    grad -= arma::vectorise(tQs * (counts - pi(current_d)) * X_pi.col(i).t());
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}
double mnhmm::objective_A(const arma::vec& x, arma::vec& grad) {
  mstep_iter++;
  double value = 0;
  arma::mat eta_Arow = arma::mat(x.memptr(), S - 1, X_A.n_rows);
  arma::mat gamma_Arow = sum_to_zero(eta_Arow, Qs);
  arma::vec A1(S);
  arma::vec log_A1(S);
  grad.zeros();
  arma::uvec idx(S);
  if (!iv_A && !tv_A) {
    A1 = softmax(gamma_Arow * X_A.slice(0).col(0));
    log_A1 = log(A1);
  }
  arma::mat tQs = Qs.t();
  for (arma::uword i = 0; i < N; ++i) {
    if (iv_A && !tv_A) {
      A1 = softmax(gamma_Arow * X_A.slice(i).col(0));
      log_A1 = log(A1);
    }
    for (arma::uword t = 1; t < Ti(i); ++t) {
      const arma::vec& counts = E_A(current_d, current_s).slice(t).col(i);
      idx = arma::find(counts);
      double sum_ea = arma::accu(counts.rows(idx));
      if (tv_A) {
        A1 = softmax(gamma_Arow * X_A.slice(i).col(t));
        log_A1 = log(A1);
      }
      double val = arma::dot(counts.rows(idx), log_A1.rows(idx));
      if (!std::isfinite(val)) {
        grad.zeros();
        return maxval;
      }
      value -= val;
      grad -= arma::vectorise(tQs * (counts - sum_ea * A1) * X_A.slice(i).col(t).t());
    }
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}
double mnhmm::objective_B(const arma::vec& x, arma::vec& grad) {
  mstep_iter++;
  double value = 0;
  arma::uword Mc = M(current_c);
  arma::mat eta_Brow = arma::mat(x.memptr(), Mc - 1, X_B(current_c).n_rows);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm(current_c));
  arma::vec B1(Mc);
  arma::vec log_B1(Mc);
  grad.zeros();
  if (!iv_B(current_c) && !tv_B(current_c)) {
    B1 = softmax(gamma_Brow * X_B(current_c).slice(0).col(0));
    log_B1 = log(B1);
  }
  arma::mat I(Mc, Mc, arma::fill::eye);
  arma::mat tQm = Qm(current_c).t();
  for (arma::uword i = 0; i < N; ++i) {
    if (iv_B(current_c) && !tv_B(current_c)) {
      B1 = softmax(gamma_Brow * X_B(current_c).slice(i).col(0));
      log_B1 = log(B1);
    }
    for (arma::uword t = 0; t < Ti(i); ++t) {
      if (obs(current_c, t, i) < Mc) {
        double e_b = E_B(current_d, current_c)(t, i, current_s);
        if (e_b > 0) {
          if (tv_B(current_c)) {
            B1 = softmax(gamma_Brow * X_B(current_c).slice(i).col(t));
            log_B1 = log(B1);
          }
          double val = e_b * log_B1(obs(current_c, t, i));
          if (!std::isfinite(val)) {
            grad.zeros();
            return maxval;
          }
          value -= val;
          grad -= arma::vectorise(tQm * e_b  * (I.col(obs(current_c, t, i)) - B1) * 
            X_B(current_c).slice(i).col(t).t());
        }
      }
    }
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}
double mnhmm::objective_omega(const arma::vec& x, arma::vec& grad) {
  
  mstep_iter++;
  double value = 0;
  
  eta_omega = arma::mat(x.memptr(), D - 1, X_omega.n_rows);
  gamma_omega = sum_to_zero(eta_omega, Qd);
  grad.zeros();
  arma::mat tQd = Qd.t();
  arma::uvec idx(D);
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    const arma::vec& counts = E_omega.col(i);
    idx = arma::find(counts);
    double val = arma::dot(counts.rows(idx), log_omega.rows(idx));
    if (!std::isfinite(val)) {
      grad.zeros();
      return maxval;
    }
    value -= val;
    grad -= arma::vectorise(tQd * (counts - omega) * X_omega.col(i).t());
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}

void mnhmm::mstep_pi(const arma::uword print_level) {
  mstep_return_code = 0;
  // Use closed form solution
  if (icpt_only_pi && lambda < 1e-12) {
    for (arma::uword d = 0; d < D; ++d) {
      eta_pi(d) = Qs.t() * log(arma::sum(E_pi(d), 1) + arma::datum::eps);
      if (!eta_pi(d).is_finite()) {
        mstep_return_code = -100;
        return;
      }
    }
    return;
  }
  if (!opt_pi) {
    Rcpp::stop("Optimizer opt_pi not initialized! Shouldn't be possible, file an issue.");
  }
  nlopt_set_min_objective(opt_pi, objective_pi_wrapper, this);
  arma::vec x(eta_pi(0).memptr(), eta_pi.n_elem, false, true);
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_pi(0).n_elem);
  mstep_iter = 0;
  for (arma::uword d = 0; d < D; ++d) {
    current_d = d;
    arma::vec x(eta_pi(d).memptr(), eta_pi(d).n_elem, false, true);
    ll = objective_pi(x, grad);
    mstep_iter = 0;
    if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
      return_code = 1; // already converged (L-BFGS gradient tolerance)
    } else {
      return_code = nlopt_optimize(opt_pi, x.memptr(), &minf);
      // nlopt_optimize can return generic failure code due to small gradients
      if (return_code == -1) {
        double ll_new = objective_pi(x, grad);
        double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
        if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < 1e-8) {
          return_code = 7;
        }
      }
    }
    if (print_level > 2) {
      Rcpp::Rcout<<"M-step of initial probabilities ended with return code "<<
        return_code<<" after "<<mstep_iter + 1<<" iterations."<<std::endl;
    }
    if (return_code < 0) {
      mstep_return_code = return_code - 110;
      return;
    }
  }
}
void mnhmm::mstep_A(const arma::uword print_level) {
  
  mstep_return_code = 0;
  // Use closed form solution
  if (icpt_only_A && lambda < 1e-12) {
    arma::vec tmp(S);
    for (arma::uword d = 0; d < D; ++d) {
      for (arma::uword s = 0; s < S; ++s) { // from
        for (arma::uword k = 0; k < S; ++k) { //to
          tmp(k) = arma::accu(E_A(d, s).row(k));
        }
        eta_A(d).slice(s).col(0) = Qs.t() * log(tmp + arma::datum::eps);
        if (!eta_A(d).slice(s).col(0).is_finite()) {
          mstep_return_code = -200;
          return;
        }
      }
    }
    return;
  }
  if (!opt_A) {
    Rcpp::stop("Optimizer opt_A not initialized! Shouldn't be possible, file an issue.");
  }
  nlopt_set_min_objective(opt_A, objective_A_wrapper, this);
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_A(0).slice(0).n_elem);
  for (arma::uword d = 0; d < D; ++d) {
    current_d = d;
    for (arma::uword s = 0; s < S; ++s) {
      current_s = s;
      arma::vec x(eta_A(d).slice(s).memptr(), eta_A(d).slice(s).n_elem, false, true);
      ll = objective_A(x, grad);
      mstep_iter = 0;
      if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
        return_code = 1; // already converged (L-BFGS gradient tolerance)
      } else {
        return_code = nlopt_optimize(opt_A, x.memptr(), &minf);
        if (return_code == -1) {
          double ll_new = objective_A(x, grad);
          double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
          if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < 1e-8) {
            return_code = 7;
          }
        }
      }
      if (print_level > 2) {
        Rcpp::Rcout<<"M-step of transition probabilities of state "<<s + 1<<
          " ended with return code "<<return_code<<" after "<<mstep_iter + 1<<
            " iterations."<<std::endl;
      }
      if (return_code < 0) {
        mstep_return_code = return_code - 210;
        return;
      }
    }
  }
}
void mnhmm::mstep_B(const arma::uword print_level) {
  mstep_return_code = 0;
  // use closed form solution
  if (arma::all(icpt_only_B) && lambda < 1e-12) {
    for (arma::uword c = 0; c < C; ++c) {
      arma::vec tmp(M(c));
      for (arma::uword d = 0; d < D; ++d) {
        for (arma::uword s = 0; s < S; ++s) {
          tmp.zeros();
          for (arma::uword i = 0; i < N; ++i) {
            for (arma::uword t = 0; t < Ti(i); ++t) {
              if (obs(c, t, i) < M(c)) {
                tmp(obs(c, t, i)) += E_B(d, c)(t, i, s);
              }
            }
          }
          eta_B(d, c).slice(s).col(0) = Qm(c).t() * log(tmp + arma::datum::eps);
          if (!eta_B(d, c).slice(s).col(0).is_finite()) {
            mstep_return_code = -300;
            return;
          }
        }
      }
    }
    return;
  }
  double minf;
  int return_code;
  double ll;
  for (arma::uword c = 0; c < C; ++c) {
    current_c = c;
    arma::vec grad(eta_B(0, c).slice(0).n_elem);
    nlopt_set_min_objective(opt_B[c], objective_B_wrapper, this);
    for (arma::uword d = 0; d < D; ++d) {
      current_d = d;
      for (arma::uword s = 0; s < S; ++s) {
        current_s = s;
        arma::vec x(eta_B(d, c).slice(s).memptr(), eta_B(d, c).slice(s).n_elem, false, true);
        ll = objective_B(x, grad);
        mstep_iter = 0;
        if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
          return_code = 1; // already converged (L-BFGS gradient tolerance)
        } else {
          return_code = nlopt_optimize(opt_B[c], x.memptr(), &minf);
          // nlopt_optimize can return generic failure code due to small gradients
          if (return_code == -1) {
            double ll_new = objective_B(x, grad);
            double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
            if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < 1e-8) {
              return_code = 7;
            }
          }
        }
        if (print_level > 2) {
          Rcpp::Rcout<<"M-step of emission probabilities of state "<<s + 1<<
            " and channel "<<c<<" ended with return code "<<
              return_code<<" after "<<mstep_iter + 1<<
                " iterations."<<std::endl;
        }
        if (return_code < 0) {
          mstep_return_code = return_code - 310;
          return;
        }
        // todo why?
        eta_B(d, c).slice(s) = arma::mat(x.memptr(), M(c) - 1, X_B(c).n_rows);
      }
    }
  }
}
void mnhmm::mstep_omega(const arma::uword print_level) {
  mstep_return_code = 0;
  // Use closed form solution
  if (icpt_only_omega && lambda < 1e-12) {
    eta_omega = Qd.t() * log(arma::sum(E_omega, 1) + arma::datum::eps);
    if (!eta_omega.is_finite()) {
      mstep_return_code = -400;
      return;
    }
    return;
  }
  if (!opt_omega) {
    Rcpp::stop("Optimizer opt_omega not initialized! Shouldn't be possible, file an issue.");
  }
  nlopt_set_min_objective(opt_omega, objective_omega_wrapper, this);
  arma::vec x(eta_omega.memptr(), eta_omega.n_elem, false, true);
  
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_omega.n_elem);
  ll = objective_omega(x, grad);
  mstep_iter = 0;
  if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
    return_code = 1; // already converged (L-BFGS gradient tolerance)
  } else {
    return_code = nlopt_optimize(opt_omega, x.memptr(), &minf);
    // nlopt_optimize can return generic failure code due to small gradients
    if (return_code == -1) {
      double ll_new = objective_omega(x, grad);
      double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
      if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < 1e-8) {
        return_code = 7;
      }
    }
  }
  if (print_level > 2) {
    Rcpp::Rcout<<"M-step of cluster probabilities ended with return code "<<
      return_code<<" after "<<mstep_iter + 1<<" iterations."<<std::endl;
  }
  if (return_code < 0) {
    mstep_return_code = return_code - 410;
  }
}

// [[Rcpp::export]]
Rcpp::List EM_LBFGS_mnhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const arma::mat& X_omega,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool icpt_only_omega,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::field<arma::mat>& eta_pi,
    const arma::field<arma::cube>& eta_A,
    const Rcpp::List& eta_B,
    const arma::mat& eta_omega,
    const double lambda,
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, const arma::uword print_level,
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m, const arma::uword print_level_m,
    const double bound) {
  mnhmm model(
      obs, Ti, M, X_pi, X_A, X_B, X_omega, 
      icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega,
      iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega, lambda
  );
  return model.EM(
    maxeval, ftol_abs, ftol_rel, xtol_abs, xtol_rel, print_level, 
    maxeval_m, ftol_abs_m, ftol_rel_m, xtol_abs_m,  xtol_rel_m, print_level_m,
    bound
  );
}

Rcpp::List mnhmm::EM(
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, const arma::uword print_level,
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m,
    const arma::uword print_level_m, const double bound) {
  
  opt_pi = nlopt_create(NLOPT_LD_LBFGS, eta_pi(0).n_elem);
  nlopt_set_xtol_abs1(opt_pi, xtol_abs);
  nlopt_set_ftol_abs(opt_pi, ftol_abs);
  nlopt_set_xtol_rel(opt_pi, xtol_rel);
  nlopt_set_ftol_rel(opt_pi, ftol_rel);
  nlopt_set_maxeval(opt_pi, maxeval);
  nlopt_set_lower_bounds1(opt_pi, -bound);
  nlopt_set_upper_bounds1(opt_pi, bound);
  nlopt_set_vector_storage(opt_pi, 10);
  
  opt_A = nlopt_create(NLOPT_LD_LBFGS, eta_A(0).slice(0).n_elem);
  nlopt_set_xtol_abs1(opt_A, xtol_abs);
  nlopt_set_ftol_abs(opt_A, ftol_abs);
  nlopt_set_xtol_rel(opt_A, xtol_rel);
  nlopt_set_ftol_rel(opt_A, ftol_rel);
  nlopt_set_maxeval(opt_A, maxeval);
  nlopt_set_lower_bounds1(opt_A, -bound);
  nlopt_set_upper_bounds1(opt_A, bound);
  nlopt_set_vector_storage(opt_A, 10);
  
  for (arma::uword c = 0; c < C; ++c) {
    opt_B[c] = nlopt_create(NLOPT_LD_LBFGS, eta_B(0, c).slice(0).n_elem);
    nlopt_set_xtol_abs1(opt_B[c], xtol_abs);
    nlopt_set_ftol_abs(opt_B[c], ftol_abs);
    nlopt_set_xtol_rel(opt_B[c], xtol_rel);
    nlopt_set_ftol_rel(opt_B[c], ftol_rel);
    nlopt_set_maxeval(opt_B[c], maxeval);
    nlopt_set_lower_bounds1(opt_B[c], -bound);
    nlopt_set_upper_bounds1(opt_B[c], bound);
    nlopt_set_vector_storage(opt_B[c], 10);
  }
  opt_omega = nlopt_create(NLOPT_LD_LBFGS, eta_omega.n_elem);
  nlopt_set_xtol_abs1(opt_omega, xtol_abs);
  nlopt_set_ftol_abs(opt_omega, ftol_abs);
  nlopt_set_xtol_rel(opt_omega, xtol_rel);
  nlopt_set_ftol_rel(opt_omega, ftol_rel);
  nlopt_set_maxeval(opt_omega, maxeval);
  nlopt_set_lower_bounds1(opt_omega, -bound);
  nlopt_set_upper_bounds1(opt_omega, bound);
  nlopt_set_vector_storage(opt_omega, 10);
  
  arma::uword n_omega = eta_omega.n_elem;
  arma::uword n_pi = eta_pi(0).n_elem;
  arma::uword n_A = eta_A(0).n_elem;
  arma::uvec n_Bc(C);
  for (arma::uword c = 0; c < C; ++c) {
    n_Bc(c) = eta_B(0, c).n_elem;
  }
  arma::uword n_B = arma::accu(n_Bc);
  arma::rowvec new_pars(n_omega + D * (n_pi + n_A + n_B));
  arma::rowvec pars(n_omega + D * (n_pi + n_A + n_B));
  
  pars.cols(0, n_omega - 1) = arma::vectorise(eta_omega).t();
  arma::uword ii = n_omega;
  for (arma::uword d = 0; d < D; ++d) {
    pars.cols(ii, ii + n_pi - 1) = arma::vectorise(eta_pi(d)).t();
    ii += n_pi;
  }
  for (arma::uword d = 0; d < D; ++d) {
    pars.cols(ii, ii + n_A - 1) = arma::vectorise(eta_A(d)).t();
    ii += n_A;
  }
  for (arma::uword c = 0; c < C; ++c) {
    for (arma::uword d = 0; d < D; ++d) {
      pars.cols(ii, ii + n_Bc(c) - 1) = arma::vectorise(eta_B(d, c)).t();
      ii += n_Bc(c);
    }
  }
  
  double relative_change = ftol_rel + 1.0;
  double absolute_change = ftol_abs + 1.0;
  double relative_x_change = xtol_rel + 1.0;
  double absolute_x_change = xtol_abs+ 1.0;
  arma::uword iter = 0;
  double ll_new;
  double ll;
  arma::mat log_alpha(S, T);
  arma::mat log_beta(S, T);
  arma::vec loglik(N);
  arma::vec loglik_i(D);
  // Initial log-likelihood
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    if (!icpt_only_omega || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    update_log_py(i);
    for (arma::uword d = 0; d < D; ++d) {
      univariate_forward(
        log_alpha, log_pi(d), log_A(d), 
        log_py.slice(d).cols(0, Ti(i) - 1)
      );
      univariate_backward(
        log_beta, log_A(d), log_py.slice(d).cols(0, Ti(i) - 1)
      );
      loglik_i(d) = logSumExp(log_alpha.col(Ti(i) - 1));
      estep_pi(i, d, log_alpha.col(0), log_beta.col(0), loglik_i(d));
      estep_A(i, d, log_alpha, log_beta, loglik_i(d));
      estep_B(i, d, log_alpha, log_beta, loglik_i(d));
    }
    loglik(i) = logSumExp(log_omega + loglik_i);
    estep_omega(i, log_omega + loglik_i, loglik(i));
  }
  ll = arma::accu(loglik);
  double penalty_term = 0.5 * lambda * std::pow(arma::norm(pars, 2), 2);
  ll -= penalty_term;
  ll /= n_obs;
  if (print_level > 0) {
    Rcpp::Rcout<<"Initial value of the log-likelihood: "<<ll<<std::endl;
    if (print_level > 1) {
      Rcpp::Rcout<<"Initial parameter values"<<std::endl;
      Rcpp::Rcout<<pars<<std::endl;
    }
  }
  // check for user interrupt every two seconds
  auto start_time = std::chrono::steady_clock::now();
  const std::chrono::seconds check_interval(2);
  
  while (relative_change > ftol_rel && absolute_change > ftol_abs &&
         absolute_x_change > xtol_abs && 
         relative_x_change > xtol_rel && iter < maxeval) {
    
    auto current_time = std::chrono::steady_clock::now();
    if (current_time - start_time >= check_interval) {
      Rcpp::checkUserInterrupt();
      start_time = current_time; // Reset the timer
    }
    
    iter++;
    ll_new = 0;
    
    mstep_pi(print_level_m);
    if (mstep_return_code != 0) {
      return mstep_error(
        mstep_return_code, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    mstep_A(print_level_m);
    if (mstep_return_code != 0) {
      return mstep_error(
        mstep_return_code, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    mstep_B(print_level_m);
    if (mstep_return_code != 0) {
      return mstep_error(
        mstep_return_code, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    mstep_omega(print_level_m);
    if (mstep_return_code != 0) {
      return mstep_error(
        mstep_return_code, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    // Update model
    update_gamma_omega();
    update_gamma_pi();
    update_gamma_A();
    update_gamma_B();
    for (arma::uword i = 0; i < N; ++i) {
      if (!icpt_only_omega || i == 0) {
        update_omega(i);
      }
      if (!icpt_only_pi || i == 0) {
        update_pi(i);
      }
      if (iv_A || i == 0) {
        update_A(i);
      }
      if (arma::any(iv_B) || i == 0) {
        update_B(i);
      }
      update_log_py(i);
      for (arma::uword d = 0; d < D; ++d) {
        univariate_forward(
          log_alpha, log_pi(d), log_A(d), 
          log_py.slice(d).cols(0, Ti(i) - 1)
        );
        univariate_backward(
          log_beta, log_A(d), log_py.slice(d).cols(0, Ti(i) - 1)
        );
        loglik_i(d) = logSumExp(log_alpha.col(Ti(i) - 1));
        estep_pi(i, d, log_alpha.col(0), log_beta.col(0), loglik_i(d));
        estep_A(i, d, log_alpha, log_beta, loglik_i(d));
        estep_B(i, d, log_alpha, log_beta, loglik_i(d));
      }
      loglik(i) = logSumExp(log_omega + loglik_i);
      estep_omega(i, log_omega + loglik_i, loglik(i));
    }
    ll_new = arma::accu(loglik);
    new_pars.cols(0, n_omega - 1) = arma::vectorise(eta_omega).t();
    ii = n_omega;
    for (arma::uword d = 0; d < D; ++d) {
      new_pars.cols(ii, ii + n_pi - 1) = arma::vectorise(eta_pi(d)).t();
      ii += n_pi;
    }
    for (arma::uword d = 0; d < D; ++d) {
      new_pars.cols(ii, ii + n_A - 1) = arma::vectorise(eta_A(d)).t();
      ii += n_A;
    }
    for (arma::uword c = 0; c < C; ++c) {
      for (arma::uword d = 0; d < D; ++d) {
        new_pars.cols(ii, ii + n_Bc(c) - 1) = arma::vectorise(eta_B(d, c)).t();
        ii += n_Bc(c);
      }
    }
    penalty_term = 0.5 * lambda * std::pow(arma::norm(new_pars, 2), 2);
    ll_new -= penalty_term;
    ll_new /= n_obs;
    relative_change = (ll_new - ll) / (std::abs(ll) + 1e-12);
    absolute_change = ll_new - ll;
    absolute_x_change = arma::max(arma::abs(new_pars - pars));
    relative_x_change = arma::norm(new_pars - pars, 1) / arma::norm(pars, 1);
    if (print_level > 0) {
      Rcpp::Rcout<<"Iteration: "<<iter<<std::endl;
      Rcpp::Rcout<<"           "<<"log-likelihood: "<<ll_new<<std::endl;
      Rcpp::Rcout<<"           "<<"relative change of log-likelihood: "<<relative_change<<std::endl;
      Rcpp::Rcout<<"           "<<"absolute change of scaled log-likelihood: "<<absolute_change<<std::endl;
      Rcpp::Rcout<<"           "<<"relative change of parameters: "<<relative_x_change<<std::endl;
      Rcpp::Rcout<<"           "<<"maximum absolute change of parameters: "<<absolute_x_change<<std::endl;
      if (print_level > 1) {
        Rcpp::Rcout << "current parameter values"<< std::endl;
        Rcpp::Rcout<<new_pars<<std::endl;
      }
    }
    ll = ll_new;
    pars = new_pars;
    if (absolute_change < -1e-8) {
      Rcpp::warning("EM algorithm encountered decreasing log-likelihood.");
    }
  }
  int return_code = 0;
  if (iter >= maxeval) {
    return_code = 5;
  } else if (relative_change < ftol_rel || absolute_change < ftol_abs) {
    return_code = 3;
  } else if (relative_x_change < xtol_rel || absolute_x_change < xtol_abs) {
    return_code = 4;
  }
  return Rcpp::List::create(
    Rcpp::Named("return_code") = return_code,
    Rcpp::Named("eta_omega") = Rcpp::wrap(eta_omega),
    Rcpp::Named("eta_pi") = Rcpp::wrap(eta_pi),
    Rcpp::Named("eta_A") = Rcpp::wrap(eta_A),
    Rcpp::Named("eta_B") = Rcpp::wrap(eta_B),
    Rcpp::Named("logLik") = ll * n_obs,
    Rcpp::Named("penalty_term") = penalty_term,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_f_change") = relative_change,
    Rcpp::Named("absolute_f_change") = absolute_change,
    Rcpp::Named("absolute_x_change") = absolute_x_change,
    Rcpp::Named("relative_x_change") = relative_x_change
  );
}
