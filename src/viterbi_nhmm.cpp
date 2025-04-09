// Viterbi algorithm for NHMMs
#include "nhmm.h"

// [[Rcpp::export]]
Rcpp::List viterbi_nhmm(
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
    const arma::field<arma::cube>& eta_B) {
  
  nhmm model(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B
  );
  arma::umat q(model.T, model.N, arma::fill::zeros);
  arma::vec logp(model.N);
  model.viterbi(q, logp);
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
