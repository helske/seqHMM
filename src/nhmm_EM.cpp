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
  
  // EM-algorithm begins
  arma::uword n_pi = eta_pi.n_elem;
  arma::uword n_A = eta_A.n_elem;
  arma::uword n_B = eta_B.n_elem;
  arma::rowvec current_pars(n_pi + n_A + n_B);
  arma::rowvec previous_pars(n_pi + n_A + n_B);
  
  previous_pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
  previous_pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
  previous_pars.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();

  double relative_change = ftol_rel + 1.0;
  double absolute_change = ftol_abs + 1.0;
  double relative_x_change = xtol_rel + 1.0;
  double absolute_x_change = xtol_abs+ 1.0;
  arma::uword iter = 0;
  double ll_new;
  double ll;
  arma::mat log_alpha(model.S, model.T);
  arma::mat log_beta(model.S, model.T);
  // Initial log-likelihood
  for (arma::uword i = 0; i < model.N; i++) {
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
    ll += ll_i;
    
    model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i);
    model.estep_A(i, log_alpha, log_beta, ll_i);
    model.estep_B(i, log_alpha, log_beta, ll_i);
  }
  while (relative_change > ftol_rel && absolute_change > ftol_abs && 
         absolute_x_change > xtol_abs && relative_x_change > xtol_rel && iter < maxeval) {
    iter++;
    ll_new = 0;
    
    // Minimize obj(E_pi, E_A, E_B, eta_pi, eta_A, eta_B, X_pi, X_A, X_B)
    // with respect to eta_pi, eta_A, eta_B
    model.mstep_pi(xtol_abs_m, ftol_abs_m, xtol_rel_m, ftol_rel_m, maxeval_m);
    model.mstep_A(xtol_abs_m, ftol_abs_m, xtol_rel_m, ftol_rel_m, maxeval_m);
    model.mstep_B(xtol_abs_m, ftol_abs_m, xtol_rel_m, ftol_rel_m, maxeval_m);
    // Update model
    model.update_gamma_pi();
    model.update_gamma_A();
    model.update_gamma_B();
    for (arma::uword i = 0; i < model.N; i++) {
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
      model.estep_pi(i, log_alpha.col(0), log_beta.col(0), ll_i);
      model.estep_A(i, log_alpha, log_beta, ll_i);
      model.estep_B(i, log_alpha, log_beta, ll_i);
    }
    
    relative_change = (ll_new - ll) / (std::abs(ll) + 1e-8);
    absolute_change = (ll_new - ll) / n_obs;
    current_pars.cols(0, n_pi - 1) = arma::vectorise(model.eta_pi).t();
    current_pars.cols(n_pi, n_pi + n_A - 1) = arma::vectorise(model.eta_A).t();
    current_pars.cols(n_pi + n_A, n_pi + n_A + n_B - 1) = arma::vectorise(model.eta_B).t();
    absolute_x_change = arma::max(arma::abs(current_pars - previous_pars));
    relative_x_change = arma::norm(current_pars - previous_pars, 1) / arma::norm(current_pars, 1);
    if (print_level > 0) {
      Rcpp::Rcout<<"Iteration: "<<iter<<std::endl;
      Rcpp::Rcout<<"           "<<"log-likelihood: "<<ll<<std::endl;
      Rcpp::Rcout<<"           "<<"relative change of log-likelihood: "<<relative_change<<std::endl;
      Rcpp::Rcout<<"           "<<"absolute change of log-likelihood: "<<absolute_change<<std::endl;
      Rcpp::Rcout<<"           "<<"relative change of parameters: "<<relative_x_change<<std::endl;
      Rcpp::Rcout<<"           "<<"maximum absolute change of parameters: "<<absolute_x_change<<std::endl;
      if (print_level > 1) {
        Rcpp::Rcout << "current parameter values"<< std::endl;
        Rcpp::Rcout<<current_pars<<std::endl;
      }
    }
    ll = ll_new;
    previous_pars = current_pars;
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
