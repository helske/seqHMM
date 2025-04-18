// log-likelihood and gradients of NHMM
#include "config.h"
#include "fanhmm.h"
#include "create_Q.h"
#include "eta_to_gamma.h"

// [[Rcpp::export]]
Rcpp::List Rcpp_log_objective_fanhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& eta_pi,
    const arma::cube& eta_A,
    const arma::field<arma::cube>& eta_B,
    const arma::vec& prior_y,
    const Rcpp::List& W_X_B) {
  
  arma::uword S = eta_A.n_slices;
  arma::uword C = obs.n_rows;
  arma::mat Qs = create_Q(S);
  arma::field<arma::mat> Qm(C);
  for (arma::uword c = 0; c < C; ++c) {
    Qm(c) = create_Q(M(c));
  }
  
  fanhmm model(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, eta_to_gamma(eta_pi, Qs), 
      eta_to_gamma(eta_A, Qs), eta_to_gamma(eta_B, Qm), prior_y, W_X_B
  );
  return model.log_objective(Qs, Qm);
}