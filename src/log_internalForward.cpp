#include "seqHMM.h"
using namespace Rcpp;

void log_internalForward(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::icube& obs, arma::cube& alpha, int threads) {
  
#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
    for (unsigned int k = 0; k < obs.n_rows; k++) {
      alpha.slice(k).col(0) = init;
      for (unsigned int r = 0; r < obs.n_slices; r++) {
        alpha.slice(k).col(0) += emission.slice(r).col(obs(k, 0, r));
      }
      for (unsigned int t = 1; t < obs.n_cols; t++) {
        for (unsigned int i = 0; i < transition.n_rows; i++) {
          alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        }
        for (int r = 0; r < obs.n_slices; r++) {
          alpha.slice(k).col(t) += emission.slice(r).col(obs(k, t, r));
        }
      }
    }
}
