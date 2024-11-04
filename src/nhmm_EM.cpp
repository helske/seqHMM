// EM algorithm for NHMMs

#include "nhmm_forward.h"
#include "nhmm_backward.h"
//#include "mstep.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
Rcpp::List EM_LBFGS_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const arma::uvec& Ti, const arma::uword n_obs,
    const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
    const double xtol_abs, const double xtol_rel, 
    const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
    const double xtol_abs_m, const double xtol_rel_m, const arma::uword print_level) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti,
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::mat log_alpha(model.S, model.T);
  arma::mat log_beta(model.S, model.T);
  
  // EM-algorithm begins
  arma::uword n_pi = eta_pi.n_elem;
  arma::uword n_A = eta_A.n_elem;
  arma::uword n_B = eta_B.n_elem;
  arma::rowvec pars(n_pi + n_A + n_B);
  if (print_level > 1) {
    pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
    pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
    pars.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();
  }
  
  arma::field<arma::vec> E_Pi(model.N);
  arma::field<arma::cube> E_A(model.N);
  arma::field<arma::cube> E_B(model.N);
  for (arma::uword i = 0; i < model.N; i++) {
    E_Pi(i) = arma::vec(model.S);
    E_A(i) = arma::cube(model.S, model.S, model.T);
    E_B(i) = arma::cube(model.S, model.M, model.T);
  }// create opt_data object(s) already here?
  double relative_change = ftol_rel + 1.0;
  double absolute_change = ftol_abs + 1.0;
  arma::uword iter = 0;
  double ll_new;
  double ll = -1e150;
  while ((relative_change > ftol_rel) && (absolute_change > ftol_abs) && (iter < maxeval)) {
    iter++;
    ll_new = 0;
    for (arma::uword i = 0; i < model.N; i++) {
      E_Pi(i).zeros();
      E_A(i).zeros();
      E_B(i).zeros();
      if (model.iv_pi || i == 0) {
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
      
      // update parameters once more even if already converged
      // Pi
      E_Pi(i) += arma::exp(log_alpha.col(0) + log_beta.col(0) - ll_i);
      // A
      for (arma::uword j = 0; j < model.S; j++) {
        for (arma::uword k = 0; k < model.S; k++) {
          for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
            E_A(i)(k, j, t) += exp(log_alpha(k, t) + model.log_A(k, j, t) + log_beta(j, t + 1) + model.log_py(j, t + 1) - ll_i);
          }
        }
      }
      
      // B
      for (arma::uword m = 0; m < model.M; m++) {
        for (arma::uword s = 0; s < model.S; s++) {
          for (arma::uword t = 0; t < model.Ti(i); t++) {
            if (m == model.obs(t, i)) {
              E_B(i)(s, m, t) += exp(log_alpha(s, t) + log_beta(s, t) - ll_i);
            }
          }
        }
      }
    }
    
    // Minimize obj(E_pi, E_A, E_B, eta_pi, eta_A, eta_B, X_pi, X_A, X_B)
    // with respect to eta_pi, eta_A, eta_B
    model.mstep_pi(
      E_Pi, xtol_abs_m, ftol_abs_m, xtol_rel_m,
      ftol_rel_m, maxeval_m
    );
    model.mstep_A(
      E_A, xtol_abs_m, ftol_abs_m, xtol_rel_m, ftol_rel_m, maxeval_m
    );
    model.mstep_B(
      E_B, xtol_abs_m, ftol_abs_m, xtol_rel_m, ftol_rel_m, maxeval_m
    );
    relative_change = (ll_new - ll) / (std::abs(ll) + 1e-8);
    absolute_change = (ll_new - ll) / n_obs;
    ll = ll_new;
    model.update_gamma_pi();
    model.update_gamma_A();
    model.update_gamma_B();
    if (print_level > 0) {
      Rcpp::Rcout<<"Iteration: "<<iter<<std::endl;
      Rcpp::Rcout<<"           "<<"log-likelihood: "<<ll<<std::endl;
      Rcpp::Rcout<<"           "<<"relative change: "<<relative_change<<std::endl;
      Rcpp::Rcout<<"           "<<"absolute change: "<<absolute_change<<std::endl;
      if (print_level > 1) {
        pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
        pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
        pars.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();
        Rcpp::Rcout << "current parameter values"<< std::endl;
        Rcpp::Rcout<<pars<<std::endl;
      }
    }
  }
  // Final log-likelihood
  ll_new = 0;
  for (arma::uword i = 0; i < model.N; i++) {
    model.update_log_py(i);
    univariate_forward_nhmm(log_alpha, model.log_Pi, model.log_A, model.log_py.cols(0, model.Ti(i) - 1));
    ll_new += logSumExp(log_alpha.col(Ti(i) - 1));
  }
  return Rcpp::List::create(
    Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
    Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
    Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
    Rcpp::Named("logLik") = ll_new,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_change") = relative_change,
    Rcpp::Named("absolute_change") = absolute_change
  );
}
