// log-likelihood and gradients of MNHMM
#include "config.h"
#include "mnhmm.h"
#include "create_Q.h"
#include "list_to_2d_field.h"
#include "eta_to_gamma.h"

// [[Rcpp::export]]
Rcpp::List Rcpp_log_objective_mnhmm(
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
    const arma::mat& eta_omega) {
  
  arma::uword S = eta_A(0).n_slices;
  arma::uword C = obs.n_rows;
  arma::uword D = eta_A.n_elem;
  arma::mat Qs = create_Q(S);
  arma::mat Qd = create_Q(D);
  arma::field<arma::mat> Qm(C);
  for (arma::uword c = 0; c < C; ++c) {
    Qm(c) = create_Q(M(c));
  }
  mnhmm model(
      obs, Ti, M, X_pi, X_A, X_B, X_omega, 
      icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega,
      iv_A, iv_B, tv_A, tv_B, 
      eta_to_gamma(eta_pi, Qs), 
      eta_to_gamma(eta_A, Qs), 
      eta_to_gamma(list_to_2d_field(eta_B), Qm, D), 
      eta_to_gamma(eta_omega, Qd)
  );
  return model.log_objective(Qs, Qm, Qd);
}
