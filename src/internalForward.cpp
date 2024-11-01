// Internal forward algorithms for HMMs and MHMMs

#include "forward_backward.h"
#include "useomp.h"

void internalForward(const arma::mat& transition_t, const arma::cube& emission, const arma::vec& init,
  const arma::ucube& obs, arma::cube& alpha, arma::mat& scales, arma::uword threads) {

#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, scales, obs, init, emission, transition_t)
    for (arma::uword k = 0; k < obs.n_slices; k++) {
      alpha.slice(k).col(0) = init;
      for (arma::uword r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) %= emission.slice(r).col(obs(r, 0, k));
      }
      scales(0, k) = 1.0 / sum(alpha.slice(k).col(0));
      alpha.slice(k).col(0) *= scales(0, k);
      for (arma::uword t = 1; t < obs.n_cols; t++) {
        alpha.slice(k).col(t) = transition_t * alpha.slice(k).col(t - 1);
        for (arma::uword r = 0; r < obs.n_rows; r++) {
          alpha.slice(k).col(t) %= emission.slice(r).col(obs(r, t, k));
        }
        scales(t, k) = 1.0 / sum(alpha.slice(k).col(t));
        alpha.slice(k).col(t) *= scales(t, k);
      }
    }
}

void internalForwardx(const arma::sp_mat& transition_t, const arma::cube& emission,
  const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, arma::mat& scales,
  arma::uword threads) {

#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, scales, obs, init, emission, transition_t)
    for (arma::uword k = 0; k < obs.n_slices; k++) {

      alpha.slice(k).col(0) = init.col(k);
      for (arma::uword r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) %= emission.slice(r).col(obs(r, 0, k));
      }
      scales(0, k) = 1.0 / sum(alpha.slice(k).col(0));
      alpha.slice(k).col(0) *= scales(0, k);
      for (arma::uword t = 1; t < obs.n_cols; t++) {
        alpha.slice(k).col(t) = transition_t * alpha.slice(k).col(t - 1);
        for (arma::uword r = 0; r < obs.n_rows; r++) {
          alpha.slice(k).col(t) %= emission.slice(r).col(obs(r, t, k));
        }
        scales(t, k) = 1.0 / sum(alpha.slice(k).col(t));
        alpha.slice(k).col(t) *= scales(t, k);
      }
    }
}
