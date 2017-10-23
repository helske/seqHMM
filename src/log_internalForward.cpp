// Internal forward algorithm for HMMs and MHMMs using log-space

#include "forward_backward.h"
#include "logsumexp.h"

void log_internalForward(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube& obs, arma::cube& alpha, unsigned int threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      alpha.slice(k).col(0) = init;
      for (unsigned int r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) += emission.slice(r).col(obs(r, 0, k));
      }
      for (unsigned int t = 1; t < obs.n_cols; t++) {
        for (unsigned int i = 0; i < transition.n_rows; i++) {
          alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        }
        for (unsigned int r = 0; r < obs.n_rows; r++) {
          alpha.slice(k).col(t) += emission.slice(r).col(obs(r, t, k));
        }
      }
    }
}

void log_internalForwardx(const arma::mat& transition, const arma::cube& emission,
                          const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, unsigned int threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      alpha.slice(k).col(0) = init.col(k);
      for (unsigned int r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) += emission.slice(r).col(obs(r, 0, k));
      }
      for (unsigned int t = 1; t < obs.n_cols; t++) {
        for (unsigned int i = 0; i < transition.n_rows; i++) {
          alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        }
        for (unsigned int r = 0; r < obs.n_rows; r++) {
          alpha.slice(k).col(t) += emission.slice(r).col(obs(r, t, k));
        }
      }
    }
}
