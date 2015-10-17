#include "seqHMM.h"
using namespace Rcpp;

void internalForwardx(const arma::mat& transition, const arma::cube& emission,
    const arma::mat& init, const arma::icube& obs, arma::cube& alpha, arma::mat& scales,
    int threads) {

#pragma omp parallel for schedule(static) num_threads(threads)
  for (unsigned int k = 0; k < obs.n_rows; k++) {
    for (unsigned int i = 0; i < emission.n_rows; i++) {
      alpha(i, 0, k) = init(i, k);
      for (unsigned int r = 0; r < obs.n_slices; r++) {
        alpha(i, 0, k) *= emission(i, obs(k, 0, r), r);
      }
    }
    scales(0, k) = sum(alpha.slice(k).col(0));
    alpha.slice(k).col(0) /= scales(0, k);
    for (unsigned int t = 1; t < obs.n_cols; t++) {
      for (unsigned int i = 0; i < transition.n_rows; i++) {
        alpha(i, t, k) = arma::dot(transition.col(i), alpha.slice(k).col(t - 1));
        for (int r = 0; r < obs.n_slices; r++) {
          alpha(i, t, k) *= emission(i, obs(k, t, r), r);
        }
      }
      scales(t, k) = sum(alpha.slice(k).col(t));
      alpha.slice(k).col(t) /= scales(t, k);
    }
  }
}
