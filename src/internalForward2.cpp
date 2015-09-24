
#include "seqHMM.h"
using namespace Rcpp;

void internalForward2(const arma::mat& transition, const arma::cube& emission, 
  const arma::vec& init, const arma::icube& obs, arma::cube& alpha) {  
  
  for(unsigned int k = 0; k < obs.n_rows; k++){      
    for(unsigned int i=0; i < emission.n_rows; i++){      
      alpha(i,0,k) = init(i);
      for(unsigned int r = 0; r < obs.n_slices; r++){
        alpha(i,0,k) *= emission(i,obs(k,0,r),r);
      }
    }
    alpha.slice(k).col(0) /= arma::accu(alpha.slice(k).col(0));
    for(unsigned int t = 1; t < obs.n_cols; t++){  
      for(unsigned int i = 0; i < transition.n_rows; i++){
        alpha(i, t, k) = arma::as_scalar(transition.col(i).t()*alpha.slice(k).col(t - 1));
        for(int r = 0; r < obs.n_slices; r++){
          alpha(i, t, k) *= emission(i,obs(k,t,r),r);
        }
      }
      alpha.slice(k).col(t) /= arma::accu(alpha.slice(k).col(t));
    }    
  }
}