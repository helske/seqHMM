#include "seqHMM.h"
using namespace Rcpp;

void log_internalForwardx(const arma::mat& transition, const arma::cube& emission,
    const arma::mat& init, const arma::icube& obs, arma::cube& alpha, int threads) {

#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
  for (unsigned int k = 0; k < obs.n_rows; k++) {
    for (unsigned int i = 0; i < emission.n_rows; i++) {
      alpha(i, 0, k) = init(i, k);
      for (unsigned int r = 0; r < obs.n_slices; r++) {
        alpha(i, 0, k) += emission(i, obs(k, 0, r), r);
      }
    }

    for (unsigned int t = 1; t < obs.n_cols; t++) {
      for (unsigned int i = 0; i < transition.n_rows; i++) {
        alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        for (int r = 0; r < obs.n_slices; r++) {
          alpha(i, t, k) += emission(i, obs(k, t, r), r);
        }
      }
    }
  }
}
