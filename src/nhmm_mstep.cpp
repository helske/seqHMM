#include "eta_to_gamma.h"
#include "get_parameters.h"
#include "logsumexp.h"
#include "sum_to_zero.h"
#include "nhmm_sc.h"
#include <nloptrAPI.h>

double nhmm_base::objective_pi(const arma::vec& x, arma::vec& grad, const nhmm_opt_data_pi& opt_data) {
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
    double sum_e = sum(opt_data.E_Pi(i));
    value -= arma::as_scalar(opt_data.E_Pi(i).t() * tmp - sum_e * logSumExp(tmp));
    grad -= Qt * (opt_data.E_Pi(i) - sum_e * Pi) * X_pi.col(i).t();
  }
  return value;
}

double nhmm_base::objective_A(const arma::vec& x, arma::vec& grad, const nhmm_opt_data_A& opt_data) {
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
      double sum_e = sum(opt_data.E_A(i).slice(t).row(opt_data.s));
      value -= arma::as_scalar(opt_data.E_A(i).slice(t).row(opt_data.s) * tmp - sum_e * logSumExp(tmp));
      grad -= arma::vectorise(Qt * (opt_data.E_A(i).slice(t).row(opt_data.s).t() - sum_e * A1) * X_A.slice(i).col(t).t());
    }
  }
  return value;
}
double nhmm_sc::objective_B(const arma::vec& x, arma::vec& grad, const nhmm_sc_opt_data_B& opt_data) {
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
    for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
      if (tv_B) {
        tmp = gamma_Brow * X_B.slice(i).col(t);
        B1 = softmax(tmp);
      }
      double sum_e = sum(opt_data.E_B(i).slice(t).row(opt_data.s));
      value -= arma::as_scalar(opt_data.E_B(i).slice(t).row(opt_data.s) * tmp - sum_e * logSumExp(tmp));
      grad -= arma::vectorise(Qt * (opt_data.E_B(i).slice(t).row(opt_data.s).t() - sum_e * B1) * X_B.slice(i).col(t).t());
    }
  }
  return value;
}
void nhmm_base::mstep_pi(const arma::field<arma::vec>& E_Pi,
                         const double xtol_abs, const double ftol_abs, const double xtol_rel,
                         const double ftol_rel, arma::uword maxeval) {
  nhmm_opt_data_pi opt_data(E_Pi);
  // Prepare wrapper data with both `this` and `opt_data` pointers
  std::pair<nhmm_base*, nhmm_opt_data_pi*> wrapper_data = {this, &opt_data};
  double minf;
  int status;
  // Pi
  // Wrap the objective function in a lambda to capture `this` pointer and `opt_data`
  auto objective_pi_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* wrapper_data = static_cast<std::pair<nhmm_sc*, nhmm_opt_data_pi*>*>(data);
    nhmm_base* self = wrapper_data->first; 
    nhmm_opt_data_pi* opt_data = wrapper_data->second;
    
    // Wrap raw pointers in Armadillo objects
    arma::vec x_vec(const_cast<double*>(x), n, false, false);
    arma::vec grad_vec(grad, n, false, true);
    
    // Call the member objective function with the provided opt_data
    return self->objective_pi(x_vec, grad_vec, *opt_data);
  };
  arma::vec x_pi = arma::vectorise(eta_pi);
  nlopt_opt opt_pi = nlopt_create(NLOPT_LD_LBFGS, x_pi.n_elem);
  nlopt_set_min_objective(opt_pi, objective_pi_wrapper, &wrapper_data);
  nlopt_set_xtol_abs1(opt_pi, xtol_abs);
  nlopt_set_ftol_abs(opt_pi, ftol_abs);
  nlopt_set_xtol_rel(opt_pi, xtol_rel);
  nlopt_set_ftol_rel(opt_pi, ftol_rel);
  nlopt_set_maxeval(opt_pi, maxeval);
  status = nlopt_optimize(opt_pi, x_pi.memptr(), &minf);
  eta_pi = arma::mat(x_pi.memptr(), S - 1, K_pi);
  nlopt_destroy(opt_pi);
}

void nhmm_base::mstep_A(const arma::field<arma::cube>& E_A, 
                         const double xtol_abs, const double ftol_abs, const double xtol_rel,
                         const double ftol_rel, arma::uword maxeval) {
  nhmm_opt_data_A opt_data(E_A);
  // Prepare wrapper data with both `this` and `opt_data` pointers
  std::pair<nhmm_base*, nhmm_opt_data_A*> wrapper_data = {this, &opt_data};
  double minf;
  int status;
  // A
  auto objective_A_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* wrapper_data = static_cast<std::pair<nhmm_base*, nhmm_opt_data_A*>*>(data);
    nhmm_base* self = wrapper_data->first; 
    nhmm_opt_data_A* opt_data = wrapper_data->second;
    // Wrap raw pointers in Armadillo objects
    arma::vec x_vec(const_cast<double*>(x), n, false, false);
    arma::vec grad_vec(grad, n, false, true);
    // Call the member objective function with the provided opt_data
    return self->objective_A(x_vec, grad_vec, *opt_data);
  };
  arma::vec x_A(eta_A.slice(0).n_elem);
  nlopt_opt opt_A = nlopt_create(NLOPT_LD_LBFGS, x_A.n_elem);
  nlopt_set_min_objective(opt_A, objective_A_wrapper, &wrapper_data);
  nlopt_set_xtol_abs1(opt_A, xtol_abs);
  nlopt_set_ftol_abs(opt_A, ftol_abs);
  nlopt_set_xtol_rel(opt_A, xtol_rel);
  nlopt_set_ftol_rel(opt_A, ftol_rel);
  nlopt_set_maxeval(opt_A, maxeval);
  for (arma::uword s = 0; s < S; s ++) {
    opt_data.s = s;
    x_A = arma::vectorise(eta_A.slice(s));
    status = nlopt_optimize(opt_A, x_A.memptr(), &minf);
    eta_A.slice(s) = arma::mat(x_A.memptr(), S - 1, K_A);
  }
  nlopt_destroy(opt_A);
}

void nhmm_sc::mstep_B(const arma::field<arma::cube>& E_B,
                         const double xtol_abs, const double ftol_abs, const double xtol_rel,
                         const double ftol_rel, arma::uword maxeval) {
  nhmm_sc_opt_data_B opt_data(E_B);
  // Prepare wrapper data with both `this` and `opt_data` pointers
  std::pair<nhmm_sc*, nhmm_sc_opt_data_B*> wrapper_data = {this, &opt_data};
  double minf;
  int status;
  // B
  auto objective_B_wrapper = [](unsigned n, const double* x, double* grad, void* data) -> double {
    auto* wrapper_data = static_cast<std::pair<nhmm_sc*, nhmm_sc_opt_data_B*>*>(data);
    nhmm_sc* self = wrapper_data->first; 
    nhmm_sc_opt_data_B* opt_data = wrapper_data->second;
    
    // Wrap raw pointers in Armadillo objects
    arma::vec x_vec(const_cast<double*>(x), n, false, false);
    arma::vec grad_vec(grad, n, false, true);
    
    // Call the member objective function with the provided opt_data
    return self->objective_B(x_vec, grad_vec, *opt_data);
  };
  arma::vec x_B(eta_B.slice(0).n_elem);
  nlopt_opt opt_B = nlopt_create(NLOPT_LD_LBFGS, x_B.n_elem);
  nlopt_set_min_objective(opt_B, objective_B_wrapper, &wrapper_data);
  nlopt_set_xtol_abs1(opt_B, xtol_abs);
  nlopt_set_ftol_abs(opt_B, ftol_abs);
  nlopt_set_xtol_rel(opt_B, xtol_rel);
  nlopt_set_ftol_rel(opt_B, ftol_rel);
  nlopt_set_maxeval(opt_B, maxeval);
  for (arma::uword s = 0; s < S; s ++) {
    opt_data.s = s;
    x_B = arma::vectorise(eta_B.slice(s));
    status = nlopt_optimize(opt_B, x_B.memptr(), &minf);
    eta_B.slice(s) = arma::mat(x_B.memptr(), M - 1, K_B);
  }
  nlopt_destroy(opt_B);
}
