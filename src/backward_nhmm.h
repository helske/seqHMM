#ifndef BACKWARD_NHMM_H
#define BACKWARD_NHMM_H

#include <RcppArmadillo.h>
#include "logsumexp.h"

template<typename submat>
void univariate_backward_nhmm(
    submat& log_beta,
    const arma::cube& log_transition, 
    const arma::mat& log_py) {
  
  unsigned int S = log_py.n_rows;
  unsigned int T = log_py.n_cols;
  
  log_beta.col(T - 1).zeros();
  for (int t = (T - 2); t >= 0; t--) {
    arma::vec tmpbeta(S);
    for (unsigned int i = 0; i < S; i++) {
      log_beta(i, t) = logSumExp(
        log_beta.col(t + 1) + log_transition.slice(t).row(i).t() + 
          log_py.col(t + 1)
      );
    }
  }
}

#endif