#include "seqHMM.h"

void internalForwardx(const arma::mat& transition, const arma::cube& emission,
    const arma::mat& init, const arma::icube& obs, arma::cube& alpha, arma::mat& scales,
    int threads) {

#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
default(none) shared(alpha, scales, obs, init, emission, transition)
  for (int k = 0; k < obs.n_rows; k++) {
    
    alpha.slice(k).col(0) = init.col(k);
    for (unsigned int r = 0; r < obs.n_slices; r++) {
      alpha.slice(k).col(0) %= emission.slice(r).col(obs(k, 0, r));
    }
    scales(0, k) = sum(alpha.slice(k).col(0));
    alpha.slice(k).col(0) /= scales(0, k);
    for (unsigned int t = 1; t < obs.n_cols; t++) {
      alpha.slice(k).col(t) = transition.t() * alpha.slice(k).col(t - 1);
      for (unsigned int r = 0; r < obs.n_slices; r++) {
        alpha.slice(k).col(t) %= emission.slice(r).col(obs(k, t, r));
      }
      scales(t, k) = sum(alpha.slice(k).col(t));
      alpha.slice(k).col(t) /= scales(t, k);
    }
  }
  if (scales.min() < 1e-100) {
    Rcpp::warning("Some of the scaling factors are smaller than 1e-100, results can be numerically unstable.");
  }
}
