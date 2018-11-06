// Internal backward algorithms for HMMs and MHMMs
#include "forward_backward.h"

void internalBackward(const arma::mat& transition, const arma::cube& emission,
    const arma::ucube& obs, arma::cube& beta, const arma::mat& scales, unsigned int threads) {

#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
default(none) shared(beta, scales, obs, emission,transition)
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    beta.slice(k).col(obs.n_cols - 1).fill(scales(obs.n_cols - 1, k));
    for (int t = obs.n_cols - 2; t >= 0; t--) {
        arma::vec tmpbeta = beta.slice(k).col(t + 1);
        for (unsigned int r = 0; r < obs.n_rows; r++) {
          tmpbeta %= emission.slice(r).col(obs(r, t + 1, k));
        }
        beta.slice(k).col(t) =  transition * tmpbeta * scales(t, k);
    }
  }
}


void internalBackwardx(const arma::sp_mat& transition, const arma::cube& emission,
  const arma::ucube& obs, arma::cube& beta, const arma::mat& scales, unsigned int threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(beta, scales, obs, emission,transition)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      beta.slice(k).col(obs.n_cols - 1).fill(scales(obs.n_cols - 1, k));
      for (int t = obs.n_cols - 2; t >= 0; t--) {
        arma::vec tmpbeta = beta.slice(k).col(t + 1);
        for (unsigned int r = 0; r < obs.n_rows; r++) {
          tmpbeta %= emission.slice(r).col(obs(r, t + 1, k));
        }
        beta.slice(k).col(t) =  transition * tmpbeta * scales(t, k);
      }
    }
}
