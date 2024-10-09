#ifndef FORWARD_NHMM_H
#define FORWARD_NHMM_H

#include <RcppArmadillo.h>
#include "logsumexp.h"

template<typename submat>
void univariate_forward_nhmm(
    submat& log_alpha,
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py) {
  
  unsigned int S = log_py.n_rows;
  unsigned int T = log_py.n_cols;
  log_alpha.col(0) = log_init + log_py.col(0);
  for (unsigned int t = 1; t < T; t++) {
    for (unsigned int i = 0; i < S; i++) {
      log_alpha(i, t) = logSumExp(
        log_alpha.col(t - 1) + log_transition.slice(t - 1).col(i) + log_py(i, t)
      );
    }
  }
}

#endif