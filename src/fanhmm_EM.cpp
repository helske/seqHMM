// EM algorithm for FAN-HMMs

#include "nhmm_forward.h"
#include "nhmm_backward.h"
#include "fanhmm_sc.h"
#include "nhmm_mc.h"
#include "logsumexp.h"
#include "sum_to_zero.h"
#include "mstep_error.h"
#include <nloptrAPI.h>
#include <chrono>

// only for case where A_t depends on y_t, otherwise call base_nhmm function
double fanhmm_sc::objective_A(const arma::vec& x, arma::vec& grad) {
  
  mstep_iter++;
  double value = 0;
  arma::mat eta_Arow = arma::mat(x.memptr(), S - 1, K_A);
  arma::mat gamma_Arow = sum_to_zero(eta_Arow, Qs);
  arma::cube rho_As = arma::cube(x.memptr() + (S - 1) * K_A, S - 1, L_A, M - 1);
  arma::uword n = (S - 1) * L_A;
  arma::cube phi_As = rho_to_phi(rho_As, Qs);
  arma::vec A1(S);
  arma::vec log_A1(S);
  grad.zeros();
  arma::mat tQs = Qs.t();
  arma::uvec idx(S);
  arma::vec tmp(S - 1);
  for (arma::uword i = 0; i < N; i++) {
    for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
      const arma::vec& counts = E_A(current_s).slice(t).col(i);
      idx = arma::find(counts);
      double sum_ea = arma::accu(counts.rows(idx));
      A1 = softmax(
        gamma_Arow * X_A.slice(i).col(t) +
          phi_As.slice(obs(t, i)) * W_A.slice(i).col(t)
      );
      log_A1 = log(A1);
      
      double val = arma::dot(counts.rows(idx), log_A1.rows(idx));
      if (!std::isfinite(val)) {
        grad.zeros();
        return maxval;
      }
      value -= val;
      tmp = tQs * (counts - sum_ea * A1);
      grad.subvec(0, (K_A * (S - 1)) - 1) -= 
        arma::vectorise(tmp * X_A.slice(i).col(t).t());
      if (obs(t, i) > 0) {
        grad.subvec(K_A * (S - 1) + (obs(t, i) - 1) * n, K_A * (S - 1) + obs(t, i) * n - 1) -= 
          arma::vectorise(tmp * W_A.slice(i).col(t).t());
      }
    }
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}
void fanhmm_sc::mstep_A(const double ftol_abs, const double ftol_rel, 
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
    auto* self = static_cast<fanhmm_sc*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_A(x_vec, grad_vec);
  };
  
  arma::uword n_eta = eta_A.slice(0).n_elem;
  arma::uword n_rho = rho_A(0).n_elem;
  arma::uword n = n_eta + n_rho;
  arma::vec x(n);
  nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, n);
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
  arma::vec grad(n);
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    x = arma::join_cols(
      arma::vectorise(eta_A.slice(s)), 
      arma::vectorise(rho_A(s))
    );
    ll = objective_A(x, grad);
    mstep_iter = 0;
    if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
      return_code = 1; // already converged (L-BFGS gradient tolerance)
    } else {
      return_code = nlopt_optimize(opt, x.memptr(), &minf);
      // nlopt_optimize can return generic failure code due to small gradients
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
    eta_A.slice(s) = arma::mat(x.memptr(), S - 1, K_A);
    rho_A(s) = arma::cube(x.memptr() + (S - 1) * K_A, S - 1, L_A, M - 1);
    if (return_code < 0) {
      mstep_return_code = return_code - 210;
      nlopt_destroy(opt);
      return;
    }
  }
  nlopt_destroy(opt);
}

double fanhmm_sc::objective_B(const arma::vec& x, arma::vec& grad) {
  mstep_iter++;
  double value = 0;
  arma::mat eta_Brow = arma::mat(x.memptr(), M - 1, K_B);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm);
  arma::cube rho_Bs = arma::cube(x.memptr() + (M - 1) * K_B, M - 1, L_B, M - 1);
  arma::cube phi_Bs = rho_to_phi(rho_Bs, Qm);
  arma::uword n = (M - 1) * L_B;
  arma::vec B1(M);
  arma::vec log_B1(M);
  grad.zeros();
  arma::mat I(M, M, arma::fill::eye);
  arma::mat tQm = Qm.t();
  arma::vec tmp(M - 1);
  for (arma::uword i = 0; i < N; i++) {
    double e_b = E_B(0, i, current_s);
    if (e_b > 0) {
      B1 = softmax(
        gamma_Brow * X_B.slice(i).col(0) +
          phi_Bs.slice(obs_0(i)) * W_B.slice(i).col(0)
      );
      log_B1 = log(B1);
      double val = e_b * log_B1(obs(0, i));
      if (!std::isfinite(val)) {
        grad.zeros();
        return maxval;
      }
      value -= val;
      tmp = tQm * e_b * (I.col(obs(0, i)) - B1);
      grad.subvec(0, (K_B * (M - 1)) - 1) -= 
        arma::vectorise(tmp * X_B.slice(i).col(0).t());
      if (obs_0(i) > 0) {
        grad.subvec(K_B * (M - 1) + (obs_0(i) - 1) * n, K_B * (M - 1) + obs_0(i) * n - 1) -= 
          arma::vectorise(tmp * W_B.slice(i).col(0).t());
      }
    }
    for (arma::uword t = 1; t < Ti(i); t++) {
      double e_b = E_B(t, i, current_s);
      if (e_b > 0) {
        B1 = softmax(
          gamma_Brow * X_B.slice(i).col(t) +
            phi_Bs.slice(obs(t - 1, i)) * W_B.slice(i).col(t)
        );
        log_B1 = log(B1);
        double val = e_b * log_B1(obs(t, i));
        if (!std::isfinite(val)) {
          grad.zeros();
          return maxval;
        }
        value -= val;
        tmp = tQm * e_b * (I.col(obs(t, i)) - B1);
        grad.subvec(0, (K_B * (M - 1)) - 1) -= 
          arma::vectorise(tmp * X_B.slice(i).col(t).t());
        if (obs(t - 1, i) > 0) {
          grad.subvec(K_B * (M - 1) + (obs(t - 1, i) - 1) * n, K_B * (M - 1) + obs(t - 1, i) * n - 1) -= 
            arma::vectorise(tmp * W_B.slice(i).col(t).t());
        }
      }
    }
  }
  grad += lambda * x;
  return value + 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
}
void fanhmm_sc::mstep_B(const double ftol_abs, const double ftol_rel, 
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
    auto* self = static_cast<fanhmm_sc*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_B(x_vec, grad_vec);
  };
  arma::uword n_eta = eta_B.slice(0).n_elem;
  arma::uword n_rho = rho_B(0).n_elem;
  arma::uword n = n_eta + n_rho;
  arma::vec x(n);
  nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, n);
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
  arma::vec grad(n);
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    x = arma::join_cols(
      arma::vectorise(eta_B.slice(s)), 
      arma::vectorise(rho_B(s))
    );
    ll = objective_B(x, grad);
    mstep_iter = 0;
    if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
      return_code = 1; // already converged (L-BFGS gradient tolerance)
    } else {
      return_code = nlopt_optimize(opt, x.memptr(), &minf);
      // nlopt_optimize can return generic failure code due to small gradients
      if (return_code == -1) {
        double ll_new = objective_B(x, grad);
        double relative_change = abs(ll_new - ll) / (std::abs(ll) + 1e-12);
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
    eta_B.slice(s) = arma::mat(x.memptr(), M - 1, K_B);
    rho_B(s) = arma::cube(x.memptr() + (M - 1) * K_B, M - 1, L_B, M - 1);
    if (return_code < 0) {
      mstep_return_code = return_code - 310;
      nlopt_destroy(opt);
      return;
    }
  }
  nlopt_destroy(opt);
}

// [[Rcpp::export]]
Rcpp::List EM_LBFGS_fanhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& X_pi,
    const arma::cube& eta_A, const arma::cube& X_A,
    const arma::cube& eta_B, const arma::cube& X_B,
    const arma::field<arma::cube>& rho_A, const arma::cube& W_A,
    const arma::field<arma::cube>& rho_B, const arma::cube& W_B,
    const arma::uvec obs_0, const arma::umat& obs, const arma::uvec& Ti,
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B, 
    const arma::uword n_obs,
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, const arma::uword print_level,
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m, 
    const arma::uword print_level_m, const double lambda, const double bound) {
  
  fanhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B, 
      obs_0, W_A, W_B, rho_A, rho_B, n_obs, lambda
  );
  
  // EM-algorithm begins
  arma::uword S = model.S;
  arma::uword n_pi = model.eta_pi.n_elem;
  arma::uword n_A = model.eta_A.n_elem;
  arma::uword n_B = model.eta_B.n_elem;
  arma::uword n_Ar = model.rho_A(0).n_elem;
  arma::uword n_Br = model.rho_B(0).n_elem;
  arma::rowvec pars_new(n_pi + n_A + n_B + S * (n_Ar + n_Br));
  arma::rowvec pars(n_pi + n_A + n_B + S * (n_Ar + n_Br));
  pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
  pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
  pars.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();
  arma::uword ii = n_pi + n_A + n_B;
  if (n_Ar > 0) {
    for (arma::uword s = 0; s < S; s++) {
      pars.cols(ii, ii + n_Ar - 1) = arma::vectorise(model.rho_A(s)).t();
      ii += n_Ar;
    }
  }
  if (n_Br > 0) {
    for (arma::uword s = 0; s < S; s++) {
      pars.cols(ii, ii + n_Br - 1) = arma::vectorise(model.rho_B(s)).t();
      ii += n_Br;
    }
  }
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
    if (model.L_A > 0) {
      model.mstep_A(
        ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
        print_level_m
      );
    } else {
      model.nhmm_base::mstep_A(
        ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
        print_level_m
      );
    }
    if (model.mstep_return_code != 0) {
      return mstep_error_fanhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    if (model.L_B > 0) {
    model.mstep_B(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
      print_level_m
    );
    } else {
      model.nhmm_sc::mstep_B(
        ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, maxeval_m, bound, 
        print_level_m
      );
    }
    if (model.mstep_return_code != 0) {
      return mstep_error_fanhmm(
        model.mstep_return_code, model, iter, relative_change, 
        absolute_change, absolute_x_change, relative_x_change);
    }
    // Update model
    
    model.update_gamma_pi();
    model.update_gamma_A();
    model.update_gamma_B();
    model.update_phi_A();
    model.update_phi_B();
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
    ii = n_pi + n_A + n_B;
    if (n_Ar > 0) {
      for (arma::uword s = 0; s < S; s++) {
        pars_new.cols(ii, ii + n_Ar - 1) = arma::vectorise(model.rho_A(s)).t();
        ii += n_Ar;
      }
    }
    if (n_Br > 0) {
      for (arma::uword s = 0; s < S; s++) {
        pars_new.cols(ii, ii + n_Br - 1) = arma::vectorise(model.rho_B(s)).t();
        ii += n_Br;
      }
    }
    
    penalty_term = 0.5 * lambda  * std::pow(arma::norm(pars_new, 2), 2);
    ll_new -= penalty_term;
    ll_new /= model.n_obs;
    
    relative_change = (ll_new - ll) / (std::abs(ll) + 1e-8);
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
    Rcpp::Named("rho_A") = Rcpp::wrap(model.rho_A),
    Rcpp::Named("rho_B") = Rcpp::wrap(model.rho_B),
    Rcpp::Named("logLik") = ll * model.n_obs,
    Rcpp::Named("penalty_term") = penalty_term,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_f_change") = relative_change,
    Rcpp::Named("absolute_f_change") = absolute_change,
    Rcpp::Named("absolute_x_change") = absolute_x_change,
    Rcpp::Named("relative_x_change") = relative_x_change
  );
}
