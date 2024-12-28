// EM algorithm for NHMMs

#include "nhmm_forward.h"
#include "nhmm_backward.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "logsumexp.h"
#include "sum_to_zero.h"
#include "mstep_error.h"
#include <nloptrAPI.h>
#include <chrono>

double nhmm_base::objective_pi(const arma::vec& x, arma::vec& grad) {
  
  mstep_iter++;
  double value = 0;
  eta_pi = arma::mat(x.memptr(), S - 1, K_pi);
  gamma_pi = sum_to_zero(eta_pi, Qs);
  arma::mat tQs = Qs.t();
  grad.zeros();
  arma::uvec idx(S);
  for (arma::uword i = 0; i < N; i++) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    const arma::vec& counts = E_pi.col(i);
    idx = arma::find(counts);
    double val = arma::dot(counts.rows(idx), log_pi.rows(idx));
    if (!std::isfinite(val)) {
      grad.zeros();
      return maxval;
    }
    value -= val;
    grad -= arma::vectorise(tQs * (counts - pi) * X_pi.col(i).t());
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}
void nhmm_base::mstep_pi(const double ftol_abs, const double ftol_rel, 
                         const double xtol_abs, const double xtol_rel, 
                         const arma::uword maxeval, const double bound, 
                         const arma::uword print_level) {
  mstep_return_code = 0;
  // Use closed form solution
  if (icpt_only_pi && lambda < 1e-12) {
    eta_pi = Qs.t() * log(arma::sum(E_pi, 1) + arma::datum::eps);
    if (!eta_pi.is_finite()) {
      mstep_return_code = -100;
      return;
    }
    return;
  }
  
  auto objective_pi_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_base*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_pi(x_vec, grad_vec);
  };
  
  arma::vec x(eta_pi.memptr(), eta_pi.n_elem, false, true);
  nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, x.n_elem);
  nlopt_set_min_objective(opt, objective_pi_wrapper, this);
  nlopt_set_xtol_abs1(opt, xtol_abs);
  nlopt_set_ftol_abs(opt, ftol_abs);
  nlopt_set_xtol_rel(opt, xtol_rel);
  nlopt_set_ftol_rel(opt, ftol_rel);
  nlopt_set_maxeval(opt, maxeval);
  nlopt_set_lower_bounds1(opt, -bound);
  nlopt_set_upper_bounds1(opt, bound);
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_pi.n_elem);
  ll = objective_pi(x, grad);
  mstep_iter = 0;
  if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
    return_code = 1; // already converged (L-BFGS gradient tolerance)
  } else {
    return_code = nlopt_optimize(opt, x.memptr(), &minf);
    // nlopt_optimize can return generic failure code due to small gradients
    if (return_code == -1) {
      double ll_new = objective_pi(x, grad);
      double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
      if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < ftol_rel) {
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
  }
  nlopt_destroy(opt);
}

double nhmm_base::objective_A(const arma::vec& x, arma::vec& grad) {
  
  mstep_iter++;
  double value = 0;
  arma::mat eta_Arow = arma::mat(x.memptr(), S - 1, K_A);
  arma::mat gamma_Arow = sum_to_zero(eta_Arow, Qs);
  arma::vec A1(S);
  arma::vec log_A1(S);
  grad.zeros();
  if (!iv_A && !tv_A) {
    A1 = softmax(gamma_Arow * X_A.slice(0).col(0));
    log_A1 = log(A1);
  }
  arma::mat tQs = Qs.t();
  arma::uvec idx(S);
  for (arma::uword i = 0; i < N; i++) {
    if (iv_A && !tv_A) {
      A1 = softmax(gamma_Arow * X_A.slice(i).col(0));
      log_A1 = log(A1);
    }
    for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
      const arma::vec& counts = E_A(current_s).slice(t).col(i);
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
void nhmm_base::mstep_A(const double ftol_abs, const double ftol_rel, 
                        const double xtol_abs, const double xtol_rel, 
                        const arma::uword maxeval, const double bound, 
                        const arma::uword print_level) {
  mstep_return_code = 0;
  
  // Use closed form solution
  if (icpt_only_A && lambda < 1e-12) {
    arma::vec tmp(S);
    for (arma::uword s = 0; s < S; s++) { // from
      for (arma::uword k = 0; k < S; k++) { //to
        tmp(k) = arma::accu(E_A(s).row(k));
      }
      eta_A.slice(s).col(0) = Qs.t() * log(tmp + arma::datum::eps);
      if (!eta_A.slice(s).col(0).is_finite()) {
        mstep_return_code = -200;
        return;
      }
    }
    return;
  }
  
  auto objective_A_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_base*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_A(x_vec, grad_vec);
  };
  
  nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, eta_A.slice(0).n_elem);
  nlopt_set_min_objective(opt, objective_A_wrapper, this);
  nlopt_set_xtol_abs1(opt, xtol_abs);
  nlopt_set_ftol_abs(opt, ftol_abs);
  nlopt_set_xtol_rel(opt, xtol_rel);
  nlopt_set_ftol_rel(opt, ftol_rel);
  nlopt_set_maxeval(opt, maxeval);
  nlopt_set_lower_bounds1(opt, -bound);
  nlopt_set_upper_bounds1(opt, bound);
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_A.slice(0).n_elem);
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    arma::vec x(eta_A.slice(s).memptr(), eta_A.slice(s).n_elem, false, true);
    ll = objective_A(x, grad);
    mstep_iter = 0;
    if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
      return_code = 1; // already converged (L-BFGS gradient tolerance)
    } else {
      return_code = nlopt_optimize(opt, x.memptr(), &minf);
      if (return_code == -1) {
        double ll_new = objective_A(x, grad);
        double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
        if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < ftol_rel) {
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
      nlopt_destroy(opt);
      return;
    }
  }
  nlopt_destroy(opt);
}

double nhmm_sc::objective_B(const arma::vec& x, arma::vec& grad) {
  mstep_iter++;
  double value = 0;
  arma::mat eta_Brow = arma::mat(x.memptr(), M - 1, K_B);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm);
  arma::vec B1(M);
  arma::vec log_B1(M);
  grad.zeros();
  if (!iv_B && !tv_B) {
    B1 = softmax(gamma_Brow * X_B.slice(0).col(0));
    log_B1 = log(B1);
  }
  arma::mat I(M, M, arma::fill::eye);
  arma::mat tQm = Qm.t();
  for (arma::uword i = 0; i < N; i++) {
    if (iv_B && !tv_B) {
      B1 = softmax(gamma_Brow * X_B.slice(i).col(0));
      log_B1 = log(B1);
    }
    for (arma::uword t = 0; t < Ti(i); t++) {
      if (obs(t, i) < M) {
        double e_b = E_B(t, i, current_s);
        if (e_b > 0) {
          if (tv_B) {
            B1 = softmax(gamma_Brow * X_B.slice(i).col(t));
            log_B1 = log(B1);
          }
          double val = e_b * log_B1(obs(t, i));
          if (!std::isfinite(val)) {
            grad.zeros();
            return maxval;
          }
          value -= val;
          grad -= arma::vectorise(tQm * e_b * (I.col(obs(t, i)) - B1) * X_B.slice(i).col(t).t());
        }
      }
    }
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}
void nhmm_sc::mstep_B(const double ftol_abs, const double ftol_rel, 
                      const double xtol_abs, const double xtol_rel, 
                      const arma::uword maxeval, const double bound,
                      const arma::uword print_level) {
  mstep_return_code = 0;
  // use closed form solution
  if (icpt_only_B && lambda < 1e-12) {
    arma::vec tmp(M);
    for (arma::uword s = 0; s < S; s++) {
      tmp.zeros();
      for (arma::uword i = 0; i < N; i++) {
        for (arma::uword t = 0; t < Ti(i); t++) {
          if (obs(t, i) < M) {
            tmp(obs(t, i)) += E_B(t, i, s);
          }
        }
      }
      eta_B.slice(s).col(0) = Qm.t() * log(tmp + arma::datum::eps);
      if (!eta_B.slice(s).col(0).is_finite()) {
        mstep_return_code = -300;
        return;
      }
    }
    return;
  }
  auto objective_B_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_sc*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_B(x_vec, grad_vec);
  };
  nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, eta_B.slice(0).n_elem);
  nlopt_set_min_objective(opt, objective_B_wrapper, this);
  nlopt_set_xtol_abs1(opt, xtol_abs);
  nlopt_set_ftol_abs(opt, ftol_abs);
  nlopt_set_xtol_rel(opt, xtol_rel);
  nlopt_set_ftol_rel(opt, ftol_rel);
  nlopt_set_maxeval(opt, maxeval);
  nlopt_set_lower_bounds1(opt, -bound);
  nlopt_set_upper_bounds1(opt, bound);
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_B.slice(0).n_elem);
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    arma::vec x(eta_B.slice(s).memptr(), eta_B.slice(s).n_elem, false, true);
    ll = objective_B(x, grad);
    mstep_iter = 0;
    if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
      return_code = 1; // already converged (L-BFGS gradient tolerance)
    } else {
      return_code = nlopt_optimize(opt, x.memptr(), &minf);
      if (return_code == -1) {
        // nlopt_optimize can return generic failure code, check if still converged
        // doesn't really matter if fails if at least one M-step is succesful (partial M-step)
        double ll_new = objective_B(x, grad);
        double relative_change = std::abs(ll_new - ll) / (std::abs(ll) + 1e-12);
        if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < ftol_rel) {
          return_code = 7;
        }
      }
    }
    if (print_level > 2) {
      Rcpp::Rcout<<"M-step of emission probabilities of state "<<s + 1<<
        " ended with return code "<<return_code<<" after "<<mstep_iter + 1<<
          " iterations."<<std::endl;
    }
    if (return_code < 0) {
      mstep_return_code = return_code - 310;
      nlopt_destroy(opt);
      return;
    }
  }
  nlopt_destroy(opt);
}

double nhmm_mc::objective_B(const arma::vec& x, arma::vec& grad) {
  mstep_iter++;
  double value = 0;
  arma::uword Mc = M(current_c);
  arma::mat eta_Brow = arma::mat(x.memptr(), Mc - 1, K_B);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm(current_c));
  arma::vec B1(Mc);
  arma::vec log_B1(Mc);
  grad.zeros();
  if (!iv_B && !tv_B) {
    B1 = softmax(gamma_Brow * X_B.slice(0).col(0));
    log_B1 = log(B1);
  }
  arma::mat I(Mc, Mc, arma::fill::eye);
  arma::mat tQm = Qm(current_c).t();
  for (arma::uword i = 0; i < N; i++) {
    if (iv_B && !tv_B) {
      B1 = softmax(gamma_Brow * X_B.slice(i).col(0));
      log_B1 = log(B1);
    }
    for (arma::uword t = 0; t < Ti(i); t++) {
      if (obs(current_c, t, i) < Mc) {
        double e_b = E_B(current_c)(t, i, current_s);
        if (e_b > 0) {
          if (tv_B) {
            B1 = softmax(gamma_Brow * X_B.slice(i).col(t));
            log_B1 = log(B1);
          }
          double val = e_b * log_B1(obs(current_c, t, i));
          if (!std::isfinite(val)) {
            grad.zeros();
            return maxval;
          }
          value -= val;
          grad -= arma::vectorise(tQm * e_b  * (I.col(obs(current_c, t, i)) - B1) * 
            X_B.slice(i).col(t).t());
        }
      }
    }
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}

void nhmm_mc::mstep_B(const double ftol_abs, const double ftol_rel, 
                      const double xtol_abs, const double xtol_rel, 
                      const arma::uword maxeval, const double bound,
                      const arma::uword print_level) {
  mstep_return_code = 0;
  // use closed form solution
  if (icpt_only_B && lambda < 1e-12) {
    for (arma::uword c = 0; c < C; c++) {
      arma::vec tmp(M(c));
      for (arma::uword s = 0; s < S; s++) {
        tmp.zeros();
        for (arma::uword i = 0; i < N; i++) {
          for (arma::uword t = 0; t < Ti(i); t++) {
            if (obs(c, t, i) < M(c)) {
              tmp(obs(c, t, i)) += E_B(c)(t, i, s);
            }
          }
        }
        eta_B(c).slice(s).col(0) = Qm(c).t() * log(tmp + arma::datum::eps);
        if (!eta_B(c).slice(s).col(0).is_finite()) {
          mstep_return_code = -300;
          return;
        }
      }
    }
    return;
  }
  auto objective_B_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_mc*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_B(x_vec, grad_vec);
  };
  double minf;
  int return_code;
  double ll;
  for (arma::uword c = 0; c < C; c++) {
    arma::vec grad(eta_B(c).slice(0).n_elem);
    nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, eta_B(c).slice(0).n_elem);
    nlopt_set_min_objective(opt, objective_B_wrapper, this);
    nlopt_set_xtol_abs1(opt, xtol_abs);
    nlopt_set_ftol_abs(opt, ftol_abs);
    nlopt_set_xtol_rel(opt, xtol_rel);
    nlopt_set_ftol_rel(opt, ftol_rel);
    nlopt_set_maxeval(opt, maxeval);
    nlopt_set_lower_bounds1(opt, -bound);
    nlopt_set_upper_bounds1(opt, bound);
    current_c = c;
    for (arma::uword s = 0; s < S; s++) {
      current_s = s;
      arma::vec x(eta_B(c).slice(s).memptr(), eta_B(c).slice(s).n_elem, false, true);
      ll = objective_B(x, grad);
      mstep_iter = 0;
      if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
        return_code = 1; // already converged (L-BFGS gradient tolerance)
      } else {
        return_code = nlopt_optimize(opt, x.memptr(), &minf);
      }
      // nlopt_optimize can return generic failure code due to small gradients
      if (return_code == -1) {
        double ll_new = objective_B(x, grad);
        double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
        if ((arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll_new)) || relative_change < ftol_rel) {
          return_code = 7;
        }
      }
      if (print_level > 2) {
        Rcpp::Rcout<<"M-step of emission probabilities of state "<<s + 1<<
          " and channel "<<c<<" ended with return code "<<return_code<<
            " after "<<mstep_iter + 1<<" iterations."<<std::endl;
      }
      if (return_code < 0) {
        mstep_return_code = return_code - 310;
        nlopt_destroy(opt);
        return;
      }
      eta_B(c).slice(s) = arma::mat(x.memptr(), M(c) - 1, K_B);
    }
    nlopt_destroy(opt);
  }
}
// [[Rcpp::export]]
Rcpp::List EM_LBFGS_nhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& x,
    const arma::cube& eta_A, const arma::cube& X_A,
    const arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec& Ti,
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B, 
    const arma::uword n_obs,
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, const arma::uword print_level,
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m, const arma::uword print_level_m,
    const double lambda, const double bound) {
  
  nhmm_sc model(
      eta_A.n_slices, x, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B, n_obs, 
      lambda
  );
  
  // EM-algorithm begins
  arma::uword n_pi = model.eta_pi.n_elem;
  arma::uword n_A = model.eta_A.n_elem;
  arma::uword n_B = model.eta_B.n_elem;
  arma::rowvec pars_new(n_pi + n_A + n_B);
  arma::rowvec pars(n_pi + n_A + n_B);
  
  pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
  pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
  pars.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();
  
  double relative_change = ftol_rel + 1.0;
  double absolute_change = ftol_abs + 1.0;
  double relative_x_change = xtol_rel + 1.0;
  double absolute_x_change = xtol_abs + 1.0;
  arma::uword iter = 0;
  double ll_new;
  double ll = 0;
  arma::mat log_alpha(model.S, model.T);
  arma::mat log_beta(model.S, model.T);
  
  // Initial log-likelihood
  for (arma::uword i = 0; i < model.N; i++) {
    if (!model.icpt_only_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    univariate_forward_nhmm(
      log_alpha, model.log_pi, model.log_A,
      model.log_py.cols(0, model.Ti(i) - 1)
    );
    univariate_backward_nhmm(
      log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
    );
    double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
    ll += ll_i;
    model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i);
    model.estep_A(i, log_alpha, log_beta, ll_i);
    model.estep_B(i, log_alpha, log_beta, ll_i);
  }
  double penalty_term = 0.5 * lambda  * std::pow(arma::norm(pars, 2), 2);
  ll -= penalty_term;
  ll /= model.n_obs;
  
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
    
    // Minimize obj(E_pi, E_A, E_B, eta_pi, eta_A, eta_B, x, X_A, X_B)
    // with respect to eta_pi, eta_A, eta_B
    model.mstep_pi(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
      print_level_m
    );
    if (model.mstep_return_code != 0) {
      return mstep_error_nhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    model.mstep_A(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
      print_level_m
    );
    if (model.mstep_return_code != 0) {
      return mstep_error_nhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    model.mstep_B(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
      print_level_m
    );
    if (model.mstep_return_code != 0) {
      return mstep_error_nhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    // Update model
    model.update_gamma_pi();
    model.update_gamma_A();
    model.update_gamma_B();
    for (arma::uword i = 0; i < model.N; i++) {
      if (!model.icpt_only_pi || i == 0) {
        model.update_pi(i);
      }
      if (model.iv_A || i == 0) {
        model.update_A(i);
      }
      if (model.iv_B || i == 0) {
        model.update_B(i);
      }
      model.update_log_py(i);
      univariate_forward_nhmm(
        log_alpha, model.log_pi, model.log_A,
        model.log_py.cols(0, model.Ti(i) - 1)
      );
      univariate_backward_nhmm(
        log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
      );
      double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
      ll_new += ll_i;
      model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i);
      model.estep_A(i, log_alpha, log_beta, ll_i);
      model.estep_B(i, log_alpha, log_beta, ll_i);
    }
    
    pars_new.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
    pars_new.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
    pars_new.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();
    
    penalty_term = 0.5 * lambda  * std::pow(arma::norm(pars_new, 2), 2);
    ll_new -= penalty_term;
    ll_new /= model.n_obs;
    
    relative_change = (ll_new - ll) / (std::abs(ll) + 1e-12);
    absolute_change = ll_new - ll;
    absolute_x_change = arma::max(arma::abs(pars_new - pars));
    relative_x_change = arma::norm(pars_new - pars, 1) / arma::norm(pars, 1);
    if (print_level > 0) {
      Rcpp::Rcout<<"Iteration: "<<iter<<std::endl;
      Rcpp::Rcout<<"           "<<"Scaled log-likelihood: "<<ll_new<<std::endl;
      Rcpp::Rcout<<"           "<<"Relative change of log-likelihood: "<<relative_change<<std::endl;
      Rcpp::Rcout<<"           "<<"Absolute change of scaled log-likelihood: "<<absolute_change<<std::endl;
      Rcpp::Rcout<<"           "<<"Relative change of parameters: "<<relative_x_change<<std::endl;
      Rcpp::Rcout<<"           "<<"Maximum absolute change of parameters: "<<absolute_x_change<<std::endl;
      if (print_level > 1) {
        Rcpp::Rcout << "Current parameter values"<< std::endl;
        Rcpp::Rcout<<pars_new<<std::endl;
      }
    }
    ll = ll_new;
    pars = pars_new;
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
    Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
    Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
    Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
    Rcpp::Named("logLik") = ll * model.n_obs,
    Rcpp::Named("penalty_term") = penalty_term,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_f_change") = relative_change,
    Rcpp::Named("absolute_f_change") = absolute_change,
    Rcpp::Named("absolute_x_change") = absolute_x_change,
    Rcpp::Named("relative_x_change") = relative_x_change
  );
}

// [[Rcpp::export]]
Rcpp::List EM_LBFGS_nhmm_multichannel(
    const arma::mat& eta_pi, const arma::mat& x,
    const arma::cube& eta_A, const arma::cube& X_A,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec& Ti,   
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B, 
    const arma::uword n_obs,
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, const arma::uword print_level,
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m, const arma::uword print_level_m,
    const double lambda, const double bound) {
  
  nhmm_mc model(
      eta_A.n_slices, x, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B, n_obs, 
      lambda
  );
  
  // EM-algorithm begins
  arma::uword n_pi = model.eta_pi.n_elem;
  arma::uword n_A = model.eta_A.n_elem;
  arma::uvec n_Bc(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    n_Bc(c) = model.eta_B(c, 0).n_elem;
  }
  arma::uword n_B = arma::accu(n_Bc);
  arma::rowvec pars_new(n_pi + n_A + n_B);
  arma::rowvec pars(n_pi + n_A + n_B);
  
  pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
  pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
  int ii = 0;
  for (arma::uword c = 0; c < model.C; c++) {
    pars.cols(n_pi + n_A + ii, n_pi + n_A + ii + n_Bc(c) - 1) = 
      arma::vectorise(model.eta_B(c)).t();
    ii += n_Bc(c);
  }
  double relative_change = ftol_rel + 1.0;
  double absolute_change = ftol_abs + 1.0;
  double relative_x_change = xtol_rel + 1.0;
  double absolute_x_change = xtol_abs+ 1.0;
  arma::uword iter = 0;
  double ll_new;
  double ll = 0;
  arma::mat log_alpha(model.S, model.T);
  arma::mat log_beta(model.S, model.T);
  
  // Initial log-likelihood
  for (arma::uword i = 0; i < model.N; i++) {
    if (!model.icpt_only_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    univariate_forward_nhmm(
      log_alpha, model.log_pi, model.log_A,
      model.log_py.cols(0, model.Ti(i) - 1)
    );
    univariate_backward_nhmm(
      log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
    );
    double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
    ll += ll_i;
    model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i);
    model.estep_A(i, log_alpha, log_beta, ll_i);
    model.estep_B(i, log_alpha, log_beta, ll_i);
  }
  double penalty_term = 0.5 * lambda  * std::pow(arma::norm(pars, 2), 2);
  ll -= penalty_term;
  ll /= model.n_obs;
  
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
    // Minimize obj(E_pi, E_A, E_B, eta_pi, eta_A, eta_B, x, X_A, X_B)
    // with respect to eta_pi, eta_A, eta_B
    model.mstep_pi(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
      print_level_m
    );
    if (model.mstep_return_code != 0) {
      return mstep_error_nhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    model.mstep_A(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
      print_level_m
    );
    if (model.mstep_return_code != 0) {
      return mstep_error_nhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    model.mstep_B(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
      print_level_m
    );
    if (model.mstep_return_code != 0) {
      return mstep_error_nhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    // Update model
    model.update_gamma_pi();
    model.update_gamma_A();
    model.update_gamma_B();
    for (arma::uword i = 0; i < model.N; i++) {
      if (!model.icpt_only_pi || i == 0) {
        model.update_pi(i);
      }
      if (model.iv_A || i == 0) {
        model.update_A(i);
      }
      if (model.iv_B || i == 0) {
        model.update_B(i);
      }
      model.update_log_py(i);
      univariate_forward_nhmm(
        log_alpha, model.log_pi, model.log_A,
        model.log_py.cols(0, model.Ti(i) - 1)
      );
      univariate_backward_nhmm(
        log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
      );
      double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
      ll_new += ll_i;
      model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i);
      model.estep_A(i, log_alpha, log_beta, ll_i);
      model.estep_B(i, log_alpha, log_beta, ll_i);
    }
    pars_new.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
    pars_new.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
    ii = 0;
    for (arma::uword c = 0; c < model.C; c++) {
      pars.cols(n_pi + n_A + ii, n_pi + n_A + ii + n_Bc(c) - 1) = 
        arma::vectorise(model.eta_B(c)).t();
      ii += n_Bc(c);
    }
    
    penalty_term = 0.5 * lambda * std::pow(arma::norm(pars_new, 2), 2);
    ll_new -= penalty_term;
    ll_new /= model.n_obs;
    
    relative_change = (ll_new - ll) / (std::abs(ll) + 1e-12);
    absolute_change = ll_new - ll;
    absolute_x_change = arma::max(arma::abs(pars_new - pars));
    relative_x_change = arma::norm(pars_new - pars, 1) / arma::norm(pars, 1);
    if (print_level > 0) {
      Rcpp::Rcout<<"Iteration: "<<iter<<std::endl;
      Rcpp::Rcout<<"           "<<"Scaled log-likelihood: "<<ll_new<<std::endl;
      Rcpp::Rcout<<"           "<<"Relative change of log-likelihood: "<<relative_change<<std::endl;
      Rcpp::Rcout<<"           "<<"Absolute change of scaled log-likelihood: "<<absolute_change<<std::endl;
      Rcpp::Rcout<<"           "<<"Relative change of parameters: "<<relative_x_change<<std::endl;
      Rcpp::Rcout<<"           "<<"Maximum absolute change of parameters: "<<absolute_x_change<<std::endl;
      if (print_level > 1) {
        Rcpp::Rcout << "Current parameter values"<< std::endl;
        Rcpp::Rcout<<pars_new<<std::endl;
      }
    }
    ll = ll_new;
    pars = pars_new;
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
    Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
    Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
    Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
    Rcpp::Named("logLik") = ll * model.n_obs,
    Rcpp::Named("penalty_term") = penalty_term,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_f_change") = relative_change,
    Rcpp::Named("absolute_f_change") = absolute_change,
    Rcpp::Named("absolute_x_change") = absolute_x_change,
    Rcpp::Named("relative_x_change") = relative_x_change
  );
}
