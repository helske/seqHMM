// state_obs_probs algorithm for NHMM
#include "state_obs_probs.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
Rcpp::List state_obs_probs_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  arma::cube state_prob(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
}

// [[Rcpp::export]]
Rcpp::List state_obs_probs_nhmm_multichannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti,
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B,
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {

  nhmm_mc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A,
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::field<arma::cube> obs_prob(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    obs_prob(c) = arma::cube(model.M(c), model.T, model.N, arma::fill::value(arma::datum::nan));
  }
  arma::cube state_prob(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
}

// [[Rcpp::export]]
Rcpp::List state_obs_probs_mnhmm_singlechannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti,
    const bool icpt_only_omega, const bool icpt_only_pi,
    const bool icpt_only_A, const bool icpt_only_B, const bool iv_A,
    const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {

  mnhmm_sc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti,
      icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B,
      iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  arma::cube state_prob(model.S * model.D, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
}

// [[Rcpp::export]]
Rcpp::List state_obs_probs_mnhmm_multichannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti,
    const bool icpt_only_omega, const bool icpt_only_pi,
    const bool icpt_only_A, const bool icpt_only_B, const bool iv_A,
    const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {
  
  mnhmm_mc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti,
      icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B,
      iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  arma::field<arma::cube> obs_prob(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    obs_prob(c) = arma::cube(model.M(c), model.T, model.N, arma::fill::value(arma::datum::nan));
  }
  arma::cube state_prob(model.S * model.D, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
}
