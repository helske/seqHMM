#ifndef BACKWARD_NHMM_H
#define BACKWARD_NHMM_H

#include "config.h"
#include "logsumexp.h"

// time-varying A
template<typename submat>
void univariate_backward(
    submat& log_beta,
    const arma::cube& log_A, 
    const arma::mat& log_py) {
  
  arma::uword S = log_py.n_rows;
  arma::uword T = log_py.n_cols;
  
  log_beta.col(T - 1).zeros();
  for (int t = (T - 2); t >= 0; t--) {
    for (arma::uword i = 0; i < S; ++i) {
      log_beta(i, t) = logSumExp(
        log_beta.col(t + 1) + log_A.slice(t + 1).row(i).t() + log_py.col(t + 1)
      );
    }
  }
}

#endif
