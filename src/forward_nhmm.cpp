// forward algorithm for NHMM
#include "nhmm.h"

// [[Rcpp::export]]
arma::cube forward_nhmm(
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
  arma::cube log_alpha(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.forward(log_alpha);
  return log_alpha;
}
