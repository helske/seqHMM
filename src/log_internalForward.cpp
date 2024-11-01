// Internal forward algorithm for HMMs and MHMMs using log-space

#include "forward_backward.h"
#include "logsumexp.h"
#include "useomp.h"

void log_internalForward(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube& obs, arma::cube& alpha, arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
    for (arma::uword k = 0; k < obs.n_slices; k++) {
      alpha.slice(k).col(0) = init;
      for (arma::uword r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) += emission.slice(r).col(obs(r, 0, k));
      }
      for (arma::uword t = 1; t < obs.n_cols; t++) {
        for (arma::uword i = 0; i < transition.n_rows; i++) {
          alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        }
        for (arma::uword r = 0; r < obs.n_rows; r++) {
          alpha.slice(k).col(t) += emission.slice(r).col(obs(r, t, k));
        }
      }
    }
}

void log_internalForwardx(const arma::mat& transition, const arma::cube& emission,
                          const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
    for (arma::uword k = 0; k < obs.n_slices; k++) {
      alpha.slice(k).col(0) = init.col(k);
      for (arma::uword r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) += emission.slice(r).col(obs(r, 0, k));
      }
      for (arma::uword t = 1; t < obs.n_cols; t++) {
        for (arma::uword i = 0; i < transition.n_rows; i++) {
          alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        }
        for (arma::uword r = 0; r < obs.n_rows; r++) {
          alpha.slice(k).col(t) += emission.slice(r).col(obs(r, t, k));
        }
      }
    }
}
