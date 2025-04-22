#ifndef FORWARD_NHMM_H
#define FORWARD_NHMM_H

#include "config.h"
#include "logsumexp.h"

template<typename submat>
void univariate_forward(
    submat& log_alpha,
    const arma::vec& log_pi, 
    const arma::cube& log_A, 
    const arma::mat& log_py) {
  
  arma::uword S = log_py.n_rows;
  arma::uword T = log_py.n_cols;
  log_alpha.col(0) = log_pi + log_py.col(0);
  for (arma::uword t = 1; t < T; ++t) {
    for (arma::uword s = 0; s < S; ++s) {
      log_alpha(s, t) = logSumExp(
        log_alpha.col(t - 1) + log_A.slice(t).col(s) + log_py(s, t)
      );
    }
  }
}
#endif
