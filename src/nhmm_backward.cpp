// backward algorithm for NHMM
#include "nhmm_backward.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
arma::cube backward_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool iv_pi, const bool iv_A, const bool iv_B, 
    const bool tv_A, const bool tv_B) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti,
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  
  arma::cube log_beta(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  backward_nhmm(model, log_beta);
  return log_beta;
}

// [[Rcpp::export]]
arma::cube backward_nhmm_multichannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti, 
    const bool iv_pi, const bool iv_A, const bool iv_B, 
    const bool tv_A, const bool tv_B) {
  
  nhmm_mc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti,
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  
  arma::cube log_beta(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  backward_nhmm(model, log_beta);
  return log_beta;
}

// [[Rcpp::export]]
arma::cube backward_mnhmm_singlechannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool iv_omega, const bool iv_pi, const bool iv_A, const bool iv_B, 
    const bool tv_A, const bool tv_B) {
  
  mnhmm_sc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti, iv_omega, 
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  
  arma::cube log_beta(
      model.S * model.D, model.T, model.N,
      arma::fill::value(arma::datum::nan)
  );
  for (arma::uword d = 0; d < model.D; d++) {
    backward_mnhmm(model, log_beta);
  }
  return log_beta;
}

// [[Rcpp::export]]
arma::cube backward_mnhmm_multichannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti, 
    const bool iv_omega, const bool iv_pi, const bool iv_A, const bool iv_B, 
    const bool tv_A, const bool tv_B) {
  
  mnhmm_mc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti, iv_omega, 
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  
  arma::cube log_beta(
      model.S * model.D, model.T, model.N, arma::fill::value(arma::datum::nan)
  );
  backward_mnhmm(model, log_beta);
  return log_beta;
}
