
#include "seqHMM.h"
using namespace Rcpp;


void internalBackward(const arma::mat& transition, const arma::cube& emission, 
  const arma::icube& obs, arma::cube& beta, const arma::mat& scales, int threads) {  
  
#pragma omp parallel for num_threads(threads)
  for(int k = 0; k < obs.n_rows; k++){
    beta.slice(k).col(obs.n_cols - 1).fill(1.0);
    for(int t = obs.n_cols - 2; t >= 0; t--){  
      for(int i = 0; i < transition.n_rows; i++){
        arma::vec tmpbeta = beta.slice(k).col(t + 1);
        for(int r = 0; r < obs.n_slices; r++){
          tmpbeta %= emission.slice(r).col(obs(k, t + 1, r));
        }
        beta(i,t,k) = arma::as_scalar(transition.row(i) * tmpbeta)/scales(t + 1,k);        
      }
    }
  }
}