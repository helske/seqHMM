#include "seqHMM.h"

void internalBackward(const arma::mat& transition, const arma::cube& emission,
    const arma::icube& obs, arma::cube& beta, const arma::mat& scales, int threads) {

#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
default(none) shared(beta, scales, obs, emission,transition)
  for (int k = 0; k < obs.n_rows; k++) {
    beta.slice(k).col(obs.n_cols - 1).fill(1.0);
    for (int t = obs.n_cols - 2; t >= 0; t--) {
        arma::vec tmpbeta = beta.slice(k).col(t + 1);
        for (unsigned int r = 0; r < obs.n_slices; r++) {
          tmpbeta %= emission.slice(r).col(obs(k, t + 1, r));
        }
        beta.slice(k).col(t) =  transition * tmpbeta / scales(t + 1, k);
    }
  }
}
