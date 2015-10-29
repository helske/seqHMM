#include "seqHMM.h"
using namespace Rcpp;

void internalForwardx2(const arma::mat& transition, const arma::cube& emission,
  const arma::mat& init, const arma::icube& obs, arma::cube& alpha, arma::mat& scales,
  int threads) {
  
#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, scales, obs, init, emission, transition)
    for (int k = 0; k < obs.n_rows; k++) {
      arma::vec tmp = init.col(k);
      for (unsigned int r = 0; r < obs.n_slices; r++) {
        tmp %= emission.slice(r).col(obs(k, 0, r));
      }
      double sumtmp = sum(tmp);
      tmp /= sumtmp;
      scales(0, k) = log(sumtmp);
      alpha.slice(k).col(0) = log(tmp) + scales(0, k);
      
      for (unsigned int t = 1; t < obs.n_cols; t++) {
        tmp = transition * tmp;
        for (unsigned int r = 0; r < obs.n_slices; r++) {
          tmp %= emission.slice(r).col(obs(k, t, r));
        }
        sumtmp = sum(tmp);
        tmp /= sumtmp;
        scales(t, k) = scales(t - 1, k) + log(sumtmp);
        alpha.slice(k).col(t) = log(tmp) + scales(t, k);
      }
    }
}
