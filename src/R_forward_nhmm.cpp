// forward algorithm for NHMM
#include "config.h"
#include "nhmm.h"
#include "list_to_field.h"

// [[Rcpp::export]]
arma::field<arma::mat> Rcpp_forward_nhmm(
    const arma::field<arma::umat>& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::field<arma::mat>& X_A,
    const Rcpp::List& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& gamma_pi,
    const arma::cube& gamma_A,
    const arma::field<arma::cube>& gamma_B) {
  
  nhmm model(
      obs, Ti, M, X_pi, X_A, matlist_to_2d_field(X_B), 
      icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, gamma_pi, gamma_A, gamma_B
  );
  return model.forward();
}
