#ifndef FORWARD_NHMM_H
#define FORWARD_NHMM_H

#include "config.h"
#include "logsumexp.h"

template<typename submat>
void univariate_forward_logspace(
    submat& log_alpha,
    const arma::vec& log_pi, 
    const arma::cube& log_A, 
    const arma::mat& log_py) {
  
  const arma::uword S = log_py.n_rows;
  const arma::uword T = log_py.n_cols;
  log_alpha.col(0) = log_pi + log_py.col(0);
  for (arma::uword t = 1; t < T; ++t) {
    for (arma::uword s = 0; s < S; ++s) {
      log_alpha(s, t) = logSumExp(
        log_alpha.col(t - 1) + log_A.slice(t).col(s) + log_py(s, t)
      );
    }
  }
}

double forward(
    arma::mat& alpha, 
    arma::vec& scales,  
    const arma::vec& pi,
     const arma::cube& A,
               const arma::mat& py,
               const arma::uword T);

template<typename submat, typename subvec>
double univariate_forward(
    submat& alpha,
    subvec& scales, 
    const arma::vec& pi,
    const arma::cube& A,
    const arma::mat& py,
    const arma::uword T) {
  
  double ll = 0;
  alpha.col(0) = pi % py.col(0);
  double x = arma::accu(alpha.col(0));
  ll += std::log(x);
  double inv_x = 1.0 / x;
  scales(0) = inv_x;
  alpha.col(0) *= inv_x;
  for (arma::uword t = 1; t < T; ++t) {
    alpha.col(t) = A.slice(t).t() * alpha.col(t - 1) % py.col(t);
    x = arma::accu(alpha.col(t));
    ll += std::log(x);
    inv_x = 1.0 / x;
    alpha.col(t) *= inv_x;
    scales(t) = inv_x;
  }
  return ll;
}
#endif
