// EM algorithm for NHMMs

#include "nhmm_forward.h"
#include "nhmm_backward.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "logsumexp.h"
#include "sum_to_zero.h"
#include <nloptrAPI.h>

static int iter;

double nhmm_base::objective_pi(const arma::vec& x, arma::vec& grad) {
  iter++;
  
  double value = 0;
  
  eta_pi = arma::mat(x.memptr(), S - 1, K_pi);
  gamma_pi = sum_to_zero(eta_pi, Qs);
  arma::mat tmpgrad(S, K_pi, arma::fill::zeros);
  
  for (arma::uword i = 0; i < N; i++) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    double val = arma::dot(E_Pi.col(i), log_Pi);
    if (!std::isfinite(val)) {
      if (!grad.is_empty()) {
        grad.fill(std::numeric_limits<double>::max());
      }
      return std::numeric_limits<double>::max();
    }
    value -= val;
    // Only update grad if it's non-empty (i.e., for gradient-based optimization)
    if (!grad.is_empty()) {
      tmpgrad -= (E_Pi.col(i) - Pi) * X_pi.col(i).t();
      if (!tmpgrad.is_finite()) {
        grad.fill(std::numeric_limits<double>::max());
        return std::numeric_limits<double>::max();
      }
    }
  }
  grad = arma::vectorise(Qs.t() * tmpgrad);
  if (!grad.is_empty()) {
    grad += lambda * x;
  }
  return value + 0.5 * lambda * arma::dot(x, x);
}
void nhmm_base::mstep_pi(const double xtol_abs, const double ftol_abs, 
                         const double xtol_rel, const double ftol_rel, 
                         const arma::uword maxeval, 
                         const arma::uword print_level) {
  
  // Use closed form solution
  if (icpt_only_pi && lambda < 1e-12) {
    eta_pi = Qs.t() * log(arma::sum(E_Pi, 1));
    if (!eta_pi.is_finite()) {
      Rcpp::stop(
        "Some of the values in gamma_pi are nonfinite likely due to zero "
        "expected initial state counts.\n" 
        "Try increasing the penalty lambda or adding pseudocounts "
        "to avoid extreme probabilities."
      );
    }
    return;
  }
  
  auto objective_pi_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_base*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    if (grad) {
      arma::vec grad_vec(grad, n, false, true);
      return self->objective_pi(x_vec, grad_vec);
    } else {
      arma::vec grad_dummy;
      return self->objective_pi(x_vec, grad_dummy);
    }
  };
  
  arma::vec x_pi = arma::vectorise(eta_pi);
  nlopt_opt opt_pi = nlopt_create(NLOPT_LD_LBFGS, x_pi.n_elem);
  //nlopt_opt opt_pi = nlopt_create(NLOPT_LN_SBPLX, x_pi.n_elem);
  nlopt_set_min_objective(opt_pi, objective_pi_wrapper, this);
  nlopt_set_xtol_abs1(opt_pi, xtol_abs);
  nlopt_set_ftol_abs(opt_pi, ftol_abs);
  nlopt_set_xtol_rel(opt_pi, xtol_rel);
  nlopt_set_ftol_rel(opt_pi, ftol_rel);
  nlopt_set_maxeval(opt_pi, maxeval);
  double minf;
  iter = 0;
  int status = nlopt_optimize(opt_pi, x_pi.memptr(), &minf);
  if (print_level > 0 && status > 0) {
    Rcpp::Rcout<<"M-step of initial probabilities ended with status "<<status<<
      " after "<<iter<<" iterations."<<std::endl;
  }
  if (status < 0) {
    Rcpp::stop(
      "M-step of initial state probabilities errored with error code %i.\n"
      "Try increasing the penalty lambda or adding pseudocounts "
      "to avoid extreme probabilities.", status
    );
  }
  eta_pi = arma::mat(x_pi.memptr(), S - 1, K_pi);
  nlopt_destroy(opt_pi);
}

double nhmm_base::objective_A(const arma::vec& x, arma::vec& grad) {
  iter++;
  double value = 0;
  
  arma::mat eta_Arow = arma::mat(x.memptr(), S - 1, K_A);
  arma::mat gamma_Arow = sum_to_zero(eta_Arow, Qs);
  arma::vec A1(S);
  arma::vec log_A1(S);
  
  arma::mat tmpgrad(S, K_A, arma::fill::zeros);
  if (!iv_A && !tv_A) {
    A1 = softmax(gamma_Arow * X_A.slice(0).col(0));
    log_A1 = log(A1);
  }
  
  for (arma::uword i = 0; i < N; i++) {
    if (iv_A && !tv_A) {
      A1 = softmax(gamma_Arow * X_A.slice(i).col(0));
      log_A1 = log(A1);
    }
    for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
      double sum_ea = sum(E_A(current_s).slice(t).col(i));
      if (sum_ea > arma::datum::eps) {
        if (tv_A) {
          A1 = softmax(gamma_Arow * X_A.slice(i).col(t));
          log_A1 = log(A1);
        }
        double val = arma::dot(E_A(current_s).slice(t).col(i), log_A1);
        if (!std::isfinite(val)) {
          if (!grad.is_empty()) {
            grad.fill(std::numeric_limits<double>::max());
          }
          return std::numeric_limits<double>::max();
        }
        value -= val;
        
        if (!grad.is_empty()) {
          tmpgrad -= sum_ea * (E_A(current_s).slice(t).col(i) / sum_ea - A1) * X_A.slice(i).col(t).t();
          if (!tmpgrad.is_finite()) {
            grad.fill(std::numeric_limits<double>::max());
            return std::numeric_limits<double>::max();
          }
        }
      }
    }
  }
  
  grad = arma::vectorise(Qs.t() * tmpgrad);
  if (!grad.is_empty()) {
    grad += lambda * x;
  }
  return value + 0.5 * lambda * arma::dot(x, x);
}
void nhmm_base::mstep_A(const double ftol_abs, const double ftol_rel, 
                        const double xtol_abs, const double xtol_rel, 
                        const arma::uword maxeval, 
                        const arma::uword print_level) {
  // Use closed form solution
  if (icpt_only_A && lambda < 1e-12) {
    arma::vec tmp(S);
    for (arma::uword s = 0; s < S; s++) { // from
      for (arma::uword k = 0; k < S; k++) { //to
        tmp(k) = arma::accu(E_A(s).row(k));
      }
      eta_A.slice(s).col(0) = Qs.t() * log(tmp);
      if (!eta_A.slice(s).col(0).is_finite()) {
        Rcpp::stop(
          "Some of the values in gamma_A are nonfinite likely due to zero "
          "expected transition counts.\n" 
          "Try increasing the penalty lambda or adding pseudocounts "
          "to avoid extreme probabilities."
        );
      }
    }
    return;
  }
  
  auto objective_A_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_base*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    if (grad) {
      arma::vec grad_vec(grad, n, false, true);
      return self->objective_A(x_vec, grad_vec);
    } else {
      arma::vec grad_dummy;
      return self->objective_A(x_vec, grad_dummy);
    }
  };
  
  arma::vec x_A(eta_A.slice(0).n_elem);
  nlopt_opt opt_A = nlopt_create(NLOPT_LD_LBFGS, x_A.n_elem);
  //nlopt_opt opt_A = nlopt_create(NLOPT_LN_SBPLX, x_A.n_elem);
  nlopt_set_min_objective(opt_A, objective_A_wrapper, this);
  nlopt_set_xtol_abs1(opt_A, xtol_abs);
  nlopt_set_ftol_abs(opt_A, ftol_abs);
  nlopt_set_xtol_rel(opt_A, xtol_rel);
  nlopt_set_ftol_rel(opt_A, ftol_rel);
  nlopt_set_maxeval(opt_A, maxeval);
  double minf;
  int status;
  
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    x_A = arma::vectorise(eta_A.slice(s));
    iter = 0;
    status = nlopt_optimize(opt_A, x_A.memptr(), &minf);
    if (print_level > 0 && status > 0) {
      Rcpp::Rcout<<"M-step of transition probabilities of state "<<s + 1<<
        " ended with status "<<status<<" after "<<iter<<
          " iterations."<<std::endl;
    }
    if (status < 0) {
      Rcpp::stop(
        "M-step of transition probabilities errored with error code %i.\n"
        "Try increasing the penalty lambda or adding pseudocounts "
        "to avoid extreme probabilities.", status
      );
    }
    eta_A.slice(s) = arma::mat(x_A.memptr(), S - 1, K_A);
  }
  nlopt_destroy(opt_A);
}

double nhmm_sc::objective_B(const arma::vec& x, arma::vec& grad) {
  iter++;
  double value = 0;
  arma::mat eta_Brow = arma::mat(x.memptr(), M - 1, K_B);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm);
  arma::vec B1(M);
  arma::vec log_B1(M);
  arma::mat tmpgrad(M, K_B, arma::fill::zeros);
  if (!iv_B && !tv_B) {
    B1 = softmax(gamma_Brow * X_B.slice(0).col(0));
    log_B1 = log(B1);
  }
  
  arma::mat I(M, M, arma::fill::eye);
  for (arma::uword i = 0; i < N; i++) {
    if (iv_B && !tv_B) {
      B1 = softmax(gamma_Brow * X_B.slice(i).col(0));
      log_B1 = log(B1);
    }
    for (arma::uword t = 0; t < Ti(i); t++) {
      if (obs(t, i) < M) {
        double e_b = E_B(t, i, current_s);
        if (e_b > arma::datum::eps) {
          if (tv_B) {
            B1 = softmax(gamma_Brow * X_B.slice(i).col(t));
            log_B1 = log(B1);
          }
          
          double val = e_b * log_B1(obs(t, i));
          if (!std::isfinite(val)) {
            if (!grad.is_empty()) {
              grad.fill(std::numeric_limits<double>::max());
            }
            return std::numeric_limits<double>::max();
          }
          value -= val;
          if (!grad.is_empty()) {
            tmpgrad -= 
              e_b * (I.col(obs(t, i)) - B1) * X_B.slice(i).col(t).t();
            if (!tmpgrad.is_finite()) {
              grad.fill(std::numeric_limits<double>::max());
              return std::numeric_limits<double>::max();
            }
          }
        }
      }
    }
  }
  grad = arma::vectorise(Qm.t() * tmpgrad);
  if (!grad.is_empty()) {
    grad += lambda * x;
  }
  return value + 0.5 * lambda * arma::dot(x, x);
}
void nhmm_sc::mstep_B(const double ftol_abs, const double ftol_rel, 
                      const double xtol_abs, const double xtol_rel, 
                      const arma::uword maxeval, 
                      const arma::uword print_level) {
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
      eta_B.slice(s).col(0) = Qm.t() * log(tmp);
      if (!eta_B.slice(s).col(0).is_finite()) {
        Rcpp::stop(
          "Some of the values in gamma_B are nonfinite likely due to zero "
          "expected emission counts.\n" 
          "Try increasing the penalty lambda or adding pseudocounts "
          "to avoid extreme probabilities."
        );
      }
    }
    return;
  }
  auto objective_B_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_sc*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    if (grad) {
      arma::vec grad_vec(grad, n, false, true);
      return self->objective_B(x_vec, grad_vec);
    } else {
      arma::vec grad_dummy;
      return self->objective_B(x_vec, grad_dummy);
    }
  };
  arma::vec x_B(eta_B.slice(0).n_elem);
  nlopt_opt opt_B = nlopt_create(NLOPT_LD_LBFGS, x_B.n_elem);
  //nlopt_opt opt_B = nlopt_create(NLOPT_LN_SBPLX, x_B.n_elem);
  nlopt_set_min_objective(opt_B, objective_B_wrapper, this);
  nlopt_set_xtol_abs1(opt_B, xtol_abs);
  nlopt_set_ftol_abs(opt_B, ftol_abs);
  nlopt_set_xtol_rel(opt_B, xtol_rel);
  nlopt_set_ftol_rel(opt_B, ftol_rel);
  nlopt_set_maxeval(opt_B, maxeval);
  double minf;
  int status;
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    x_B = arma::vectorise(eta_B.slice(s));
    iter = 0;
    status = nlopt_optimize(opt_B, x_B.memptr(), &minf);
    if (print_level > 0 && status > 0) {
      Rcpp::Rcout<<"M-step of emission probabilities of state "<<s + 1<<
        " ended with status "<<status<<" after "<<iter<<
          " iterations."<<std::endl;
    }
    if (status < 0) {
      Rcpp::stop(
        "M-step of emission probabilities errored with error code %i.\n"
        "Try increasing the penalty lambda or adding pseudocounts "
        "to avoid extreme probabilities.", status
      );
    }
    eta_B.slice(s) = arma::mat(x_B.memptr(), M - 1, K_B);
  }
  nlopt_destroy(opt_B);
}

double nhmm_mc::objective_B(const arma::vec& x, arma::vec& grad) {
  iter++;
  double value = 0;
  arma::uword Mc = M(current_c);
  arma::mat eta_Brow = arma::mat(x.memptr(), Mc - 1, K_B);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm(current_c));
  arma::vec B1(Mc);
  arma::vec log_B1(Mc);
  arma::mat tmpgrad(Mc, K_B, arma::fill::zeros);
  
  if (!iv_B && !tv_B) {
    B1 = softmax(gamma_Brow * X_B.slice(0).col(0));
    log_B1 = log(B1);
  }
  arma::mat I(Mc, Mc, arma::fill::eye);
  
  for (arma::uword i = 0; i < N; i++) {
    if (iv_B && !tv_B) {
      B1 = softmax(gamma_Brow * X_B.slice(i).col(0));
      log_B1 = log(B1);
    }
    for (arma::uword t = 0; t < Ti(i); t++) {
      
      if (obs(current_c, t, i) < Mc) {
        double e_b = E_B(current_c)(t, i, current_s);
        if (e_b > arma::datum::eps) {
          if (tv_B) {
            B1 = softmax(gamma_Brow * X_B.slice(i).col(t));
            log_B1 = log(B1);
          }
          double val = e_b * log_B1(obs(current_c, t, i));
          if (!std::isfinite(val)) {
            if (!grad.is_empty()) {
              grad.fill(std::numeric_limits<double>::max());
            }
            return std::numeric_limits<double>::max();
          }
          value -= val;
          if (!grad.is_empty()) {
            tmpgrad -= e_b  * (I.col(obs(current_c, t, i)) - B1) * 
              X_B.slice(i).col(t).t();
            if (!tmpgrad.is_finite()) {
              grad.fill(std::numeric_limits<double>::max());
              return std::numeric_limits<double>::max();
            }
          }
        }
      }
    }
  }
  grad = arma::vectorise(Qm(current_c).t() * tmpgrad);
  if (!grad.is_empty()) {
    grad += lambda * x;
  }
  return value + 0.5 * lambda * arma::dot(x, x);
}

void nhmm_mc::mstep_B(const double ftol_abs, const double ftol_rel, 
                      const double xtol_abs, const double xtol_rel, 
                      const arma::uword maxeval, 
                      const arma::uword print_level) {
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
        eta_B(c).slice(s).col(0) = Qm(c).t() * log(tmp);
        if (!eta_B(c).slice(s).col(0).is_finite()) {
          Rcpp::stop(
            "Some of the values in gamma_B are nonfinite likely due to zero "
            "expected emission counts.\n" 
            "Try increasing the penalty lambda or adding pseudocounts "
            "to avoid extreme probabilities."
          );
        }
      }
    }
    return;
  }
  auto objective_B_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_mc*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    if (grad) {
      arma::vec grad_vec(grad, n, false, true);
      return self->objective_B(x_vec, grad_vec);
    } else {
      arma::vec grad_dummy;
      return self->objective_B(x_vec, grad_dummy);
    }
  };
  double minf;
  int status;
  for (arma::uword c = 0; c < C; c++) {
    arma::vec x_B(eta_B(c).slice(0).n_elem);
    nlopt_opt opt_B = nlopt_create(NLOPT_LD_LBFGS, x_B.n_elem);
    nlopt_set_min_objective(opt_B, objective_B_wrapper, this);
    nlopt_set_xtol_abs1(opt_B, xtol_abs);
    nlopt_set_ftol_abs(opt_B, ftol_abs);
    nlopt_set_xtol_rel(opt_B, xtol_rel);
    nlopt_set_ftol_rel(opt_B, ftol_rel);
    nlopt_set_maxeval(opt_B, maxeval);
    current_c = c;
    for (arma::uword s = 0; s < S; s++) {
      current_s = s;
      x_B = arma::vectorise(eta_B(c).slice(s));
      iter = 0;
      status = nlopt_optimize(opt_B, x_B.memptr(), &minf);
      if (print_level > 0 && status > 0) {
        Rcpp::Rcout<<"M-step of emission probabilities of state "<<s + 1<<
          " and channel "<<c<<" ended with status "<<status<<" after "<<iter<<
            " iterations."<<std::endl;
      }
      if (status < 0) {
        Rcpp::stop(
          "M-step of emission probabilities errored with error code %i.\n"
          "Try increasing the penalty lambda or adding pseudocounts "
          "to avoid extreme probabilities.", status
        );
      }
      eta_B(c).slice(s) = arma::mat(x_B.memptr(), M(c) - 1, K_B);
    }
    nlopt_destroy(opt_B);
  }
}
// [[Rcpp::export]]
Rcpp::List EM_LBFGS_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec& Ti,
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B, 
    const arma::uword n_obs,
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, const arma::uword print_level,
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m, const arma::uword print_level_m,
    const double lambda, const double pseudocount) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B, lambda
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
      log_alpha, model.log_Pi, model.log_A,
      model.log_py.cols(0, model.Ti(i) - 1)
    );
    
    univariate_backward_nhmm(
      log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
    );
    
    double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
    ll += ll_i;
    model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i, pseudocount);
    model.estep_A(i, log_alpha, log_beta, ll_i, pseudocount);
    model.estep_B(i, log_alpha, log_beta, ll_i, pseudocount);
  }
  double penalty_term = 0.5 * lambda  * arma::dot(pars, pars);
  ll -= penalty_term;
  ll /= model.n_obs;
  
  if (print_level > 0) {
    Rcpp::Rcout<<"Initial value of the log-likelihood: "<<ll<<std::endl;
    if (print_level > 1) {
      Rcpp::Rcout<<"Initial parameter values"<<std::endl;
      Rcpp::Rcout<<pars<<std::endl;
    }
  }
  while (relative_change > ftol_rel && absolute_change > ftol_abs && 
         absolute_x_change > xtol_abs && relative_x_change > xtol_rel && iter < maxeval) {
    iter++;
    ll_new = 0;
    
    // Minimize obj(E_pi, E_A, E_B, eta_pi, eta_A, eta_B, X_pi, X_A, X_B)
    // with respect to eta_pi, eta_A, eta_B
    model.mstep_pi(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_abs_m, maxeval_m, print_level_m
    );
    model.mstep_A(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_abs_m, maxeval_m, print_level_m
    );
    model.mstep_B(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_abs_m, maxeval_m, print_level_m
    );
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
        log_alpha, model.log_Pi, model.log_A,
        model.log_py.cols(0, model.Ti(i) - 1)
      );
      univariate_backward_nhmm(
        log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
      );
      double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
      ll_new += ll_i;
      model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i, pseudocount);
      model.estep_A(i, log_alpha, log_beta, ll_i, pseudocount);
      model.estep_B(i, log_alpha, log_beta, ll_i, pseudocount);
    }
    
    pars_new.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
    pars_new.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
    pars_new.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();
    
    penalty_term = 0.5 * lambda  * arma::dot(pars_new, pars_new);
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
    if (absolute_change < -1e6) {
      Rcpp::warning("EM algorithm encountered decreasing log-likelihood.");
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
    Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
    Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
    Rcpp::Named("penalized_logLik") = ll * model.n_obs,
    Rcpp::Named("penalty_term") = penalty_term * model.n_obs,
    Rcpp::Named("logLik") = (ll + penalty_term) * model.n_obs,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_f_change") = relative_change,
    Rcpp::Named("absolute_f_change") = absolute_change,
    Rcpp::Named("absolute_x_change") = absolute_x_change,
    Rcpp::Named("relative_x_change") = relative_x_change
  );
}

// [[Rcpp::export]]
Rcpp::List EM_LBFGS_nhmm_multichannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec& Ti,   
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B, 
    const arma::uword n_obs,
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, const arma::uword print_level,
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m, const arma::uword print_level_m,
    const double lambda, const double pseudocount) {
  
  nhmm_mc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B, lambda
  );
  
  // EM-algorithm begins
  arma::uword n_pi = model.eta_pi.n_elem;
  arma::uword n_A = model.eta_A.n_elem;
  arma::uvec n_Bc(model.C + 1);
  for (arma::uword c = 1; c < model.C + 1; c++) {
    n_Bc(c) = n_Bc(c - 1) + model.eta_B(c - 1).n_elem;
  }
  arma::uword n_B = n_Bc(model.C);
  arma::rowvec pars_new(n_pi + n_A + n_B);
  arma::rowvec pars(n_pi + n_A + n_B);
  
  pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
  pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
  for (arma::uword c = 0; c < model.C; c++) {
    pars.cols(
      n_pi + n_A + n_Bc(c), n_pi + n_A + n_Bc(c + 1) - 1
    ) = arma::vectorise(model.eta_B(c)).t();
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
      log_alpha, model.log_Pi, model.log_A,
      model.log_py.cols(0, model.Ti(i) - 1)
    );
    univariate_backward_nhmm(
      log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
    );
    double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
    ll += ll_i;
    model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i, pseudocount);
    model.estep_A(i, log_alpha, log_beta, ll_i, pseudocount);
    model.estep_B(i, log_alpha, log_beta, ll_i, pseudocount);
  }
  double penalty_term = 0.5 * lambda  * arma::dot(pars, pars);
  ll -= penalty_term;
  ll /= model.n_obs;
  
  if (print_level > 0) {
    Rcpp::Rcout<<"Initial value of the log-likelihood: "<<ll<<std::endl;
    if (print_level > 1) {
      Rcpp::Rcout<<"Initial parameter values"<<std::endl;
      Rcpp::Rcout<<pars<<std::endl;
    }
  }
  while (relative_change > ftol_rel && absolute_change > ftol_abs && 
         absolute_x_change > xtol_abs && relative_x_change > xtol_rel && iter < maxeval) {
    iter++;
    ll_new = 0;
    // Minimize obj(E_pi, E_A, E_B, eta_pi, eta_A, eta_B, X_pi, X_A, X_B)
    // with respect to eta_pi, eta_A, eta_B
    model.mstep_pi(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_abs_m, maxeval_m, print_level_m
    );
    model.mstep_A(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_abs_m, maxeval_m, print_level_m
    );
    model.mstep_B(
      ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_abs_m, maxeval_m, print_level_m
    );
    
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
        log_alpha, model.log_Pi, model.log_A,
        model.log_py.cols(0, model.Ti(i) - 1)
      );
      univariate_backward_nhmm(
        log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
      );
      double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
      ll_new += ll_i;
      model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i, pseudocount);
      model.estep_A(i, log_alpha, log_beta, ll_i, pseudocount);
      model.estep_B(i, log_alpha, log_beta, ll_i, pseudocount);
    }
    pars_new.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
    pars_new.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
    for (arma::uword c = 0; c < model.C; c++) {
      pars_new.cols(
        n_pi + n_A + n_Bc(c), n_pi + n_A + n_Bc(c + 1) - 1
      ) = arma::vectorise(model.eta_B(c)).t();
    }
    
    penalty_term = 0.5 * lambda * arma::dot(pars_new, pars_new);
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
    if (absolute_change < -1e6) {
      Rcpp::warning("EM algorithm encountered decreasing log-likelihood.");
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
    Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
    Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
    Rcpp::Named("penalized_logLik") = ll * model.n_obs,
    Rcpp::Named("penalty_term") = penalty_term * model.n_obs,
    Rcpp::Named("logLik") = (ll + penalty_term) * model.n_obs,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_f_change") = relative_change,
    Rcpp::Named("absolute_f_change") = absolute_change,
    Rcpp::Named("absolute_x_change") = absolute_x_change,
    Rcpp::Named("relative_x_change") = relative_x_change
  );
}
