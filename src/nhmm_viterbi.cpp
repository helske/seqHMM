// Viterbi algorithm for NHMMs
#include "nhmm_viterbi.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
Rcpp::List viterbi_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::umat q(model.T, model.N, arma::fill::zeros);
  arma::vec logp(model.N);
  viterbi_nhmm(model, q, logp);
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
// [[Rcpp::export]]
Rcpp::List viterbi_nhmm_multichannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B) {
  
  nhmm_mc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti,
      icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::umat q(model.T, model.N, arma::fill::zeros);
  arma::vec logp(model.N);
  viterbi_nhmm(model, q, logp);
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
// [[Rcpp::export]]
Rcpp::List viterbi_mnhmm_singlechannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, const bool icpt_only_omega, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B) {
  
  mnhmm_sc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti, 
      icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, 
      tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  arma::umat q(model.T, model.N, arma::fill::zeros);
  arma::vec logp(model.N);
  viterbi_mnhmm(model, q, logp);
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
// [[Rcpp::export]]
Rcpp::List viterbi_mnhmm_multichannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti, const bool icpt_only_omega, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B) {
  
  mnhmm_mc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti, 
      icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, 
      tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  arma::umat q(model.T, model.N, arma::fill::zeros);
  arma::vec logp(model.N);
  viterbi_mnhmm(model, q, logp);
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
