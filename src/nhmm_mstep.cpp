#include "eta_to_gamma.h"
#include "get_parameters.h"
#include "logsumexp.h"
#include "sum_to_zero.h"
#include "nhmm_sc.h"
#include <nloptrAPI.h>

double nhmm_base::objective_pi(const arma::vec& x, arma::vec& grad) {
  eta_pi = arma::mat(x.memptr(), S - 1, K_pi);
  gamma_pi = sum_to_zero(eta_pi, Qs);
  double value = 0;
  arma::mat Qt = Qs.t();
  grad.zeros();
  arma::vec tmp(S);
  for (arma::uword i = 0; i < N; i++) {
    if (iv_pi || i == 0) {
      tmp = gamma_pi * X_pi.col(i);
      Pi = softmax(tmp);
    }
    double sum_e = sum(E_Pi.col(i));
    value -= arma::as_scalar(E_Pi.col(i).t() * tmp - sum_e * logSumExp(tmp));
    grad -= Qt * (E_Pi.col(i) - sum_e * Pi) * X_pi.col(i).t();
  }
  return value;
}

double nhmm_base::objective_A(const arma::vec& x, arma::vec& grad) {
  arma::mat eta_Arow = arma::mat(x.memptr(), S - 1, K_A);
  arma::mat gamma_Arow = sum_to_zero(eta_Arow, Qs);
  arma::vec A1(S);
  double value = 0;
  arma::mat Qt = Qs.t();
  grad.zeros();
  arma::vec tmp(S);
  if (!iv_A && !tv_A) {
    tmp = gamma_Arow * X_A.slice(0).col(0);
    A1 = softmax(tmp);
  }
  for (arma::uword i = 0; i < N; i++) {
    if (iv_A && !tv_A) {
      tmp = gamma_Arow * X_A.slice(i).col(0);
      A1 = softmax(tmp);
    }
    for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
      if (tv_A) {
        tmp = gamma_Arow * X_A.slice(i).col(t);
        A1 = softmax(tmp);
      }
      double sum_e = sum(E_A(current_s).slice(t).col(i));
      value -= arma::as_scalar(E_A(current_s).slice(t).col(i).t() * tmp - sum_e * logSumExp(tmp));
      grad -= arma::vectorise(Qt * (E_A(current_s).slice(t).col(i) - sum_e * A1) * X_A.slice(i).col(t).t());
    }
  }
  return value;
}
double nhmm_sc::objective_B(const arma::vec& x, arma::vec& grad) {
  arma::mat eta_Brow = arma::mat(x.memptr(), M - 1, K_B);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm);
  arma::vec B1(M);
  double value = 0;
  arma::mat Qt = Qm.t();
  grad.zeros();
  arma::vec tmp(M);
  if (!iv_B && !tv_B) {
    tmp = gamma_Brow * X_B.slice(0).col(0);
    B1 = softmax(tmp);
  }
  for (arma::uword i = 0; i < N; i++) {
    if (iv_B && !tv_B) {
      tmp = gamma_Brow * X_B.slice(i).col(0);
      B1 = softmax(tmp);
    }
    for (arma::uword t = 0; t < Ti(i); t++) {
      if (tv_B) {
        tmp = gamma_Brow * X_B.slice(i).col(t);
        B1 = softmax(tmp);
      }
      double sum_e = sum(E_B(current_s).slice(t).col(i));
      value -= arma::as_scalar(E_B(current_s).slice(t).col(i).t() * tmp - sum_e * logSumExp(tmp));
      grad -= arma::vectorise(Qt * (E_B(current_s).slice(t).col(i) - sum_e * B1) * X_B.slice(i).col(t).t());
    }
  }
  return value;
}
void nhmm_base::mstep_pi(const double xtol_abs, const double ftol_abs, 
                         const double xtol_rel,
                         const double ftol_rel, arma::uword maxeval) {
  double minf;
  int status;
  // Wrap the objective function in a lambda to capture `this` pointer
  auto objective_pi_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_base*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, false);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_pi(x_vec, grad_vec);
  };
  
  arma::vec x_pi = arma::vectorise(eta_pi);
  nlopt_opt opt_pi = nlopt_create(NLOPT_LD_LBFGS, x_pi.n_elem);
  nlopt_set_min_objective(opt_pi, objective_pi_wrapper, this);
  nlopt_set_xtol_abs1(opt_pi, xtol_abs);
  nlopt_set_ftol_abs(opt_pi, ftol_abs);
  nlopt_set_xtol_rel(opt_pi, xtol_rel);
  nlopt_set_ftol_rel(opt_pi, ftol_rel);
  nlopt_set_maxeval(opt_pi, maxeval);
  status = nlopt_optimize(opt_pi, x_pi.memptr(), &minf);
  if (status < 0) {
    Rcpp::stop("M-step of initial probabilities errored with error code %i.", status);
  }
  eta_pi = arma::mat(x_pi.memptr(), S - 1, K_pi);
  nlopt_destroy(opt_pi);
}

void nhmm_base::mstep_A(const double xtol_abs, const double ftol_abs, 
                        const double xtol_rel,
                         const double ftol_rel, arma::uword maxeval) {
 
  double minf;
  int status;
  auto objective_A_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_base*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, false);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_A(x_vec, grad_vec);
  };
  
  arma::vec x_A(eta_A.slice(0).n_elem);
  nlopt_opt opt_A = nlopt_create(NLOPT_LD_LBFGS, x_A.n_elem);
  nlopt_set_min_objective(opt_A, objective_A_wrapper, this);
  nlopt_set_xtol_abs1(opt_A, xtol_abs);
  nlopt_set_ftol_abs(opt_A, ftol_abs);
  nlopt_set_xtol_rel(opt_A, xtol_rel);
  nlopt_set_ftol_rel(opt_A, ftol_rel);
  nlopt_set_maxeval(opt_A, maxeval);
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    x_A = arma::vectorise(eta_A.slice(s));
    status = nlopt_optimize(opt_A, x_A.memptr(), &minf);
    if (status < 0) {
      Rcpp::stop("M-step of transition probabilities errored with error code %i.", status);
    }
    eta_A.slice(s) = arma::mat(x_A.memptr(), S - 1, K_A);
  }
  nlopt_destroy(opt_A);
}

void nhmm_sc::mstep_B(const double xtol_abs, const double ftol_abs, 
                      const double xtol_rel,
                         const double ftol_rel, arma::uword maxeval) {
  double minf;
  int status;
  auto objective_B_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* self = static_cast<nhmm_sc*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, false);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_B(x_vec, grad_vec);
  };
  arma::vec x_B(eta_B.slice(0).n_elem);
  nlopt_opt opt_B = nlopt_create(NLOPT_LD_LBFGS, x_B.n_elem);
  nlopt_set_min_objective(opt_B, objective_B_wrapper, this);
  nlopt_set_xtol_abs1(opt_B, xtol_abs);
  nlopt_set_ftol_abs(opt_B, ftol_abs);
  nlopt_set_xtol_rel(opt_B, xtol_rel);
  nlopt_set_ftol_rel(opt_B, ftol_rel);
  nlopt_set_maxeval(opt_B, maxeval);
  for (arma::uword s = 0; s < S; s++) {
    current_s = s;
    x_B = arma::vectorise(eta_B.slice(s));
    status = nlopt_optimize(opt_B, x_B.memptr(), &minf);
    if (status < 0) {
      Rcpp::stop("M-step of emission probabilities errored with error code %i.", status);
    }
    eta_B.slice(s) = arma::mat(x_B.memptr(), M - 1, K_B);
  }
  nlopt_destroy(opt_B);
}
