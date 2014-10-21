
#include "seqHMM.h"
using namespace Rcpp;



void internalBackwardMC(const arma::mat& transition, const arma::cube& emission, 
                       const arma::Cube<int>& obs, arma::cube& beta) {  
  
  double sumtmp;
  double tmp;
  double neginf = -arma::math::inf();
  for(unsigned int k = 0; k < obs.n_rows; k++){
    beta.slice(k).col(beta.n_cols-1).fill(0.0);
    for(int t = (obs.n_cols-2); t >= 0; t--){  
      for(unsigned int i = 0; i < beta.n_rows; i++){
        sumtmp = neginf;
        for(unsigned int j = 0; j < beta.n_rows; j++){
          tmp = beta(j,t+1,k) + transition(i,j);
          for(unsigned int r = 0; r < obs.n_slices; r++){
            tmp += emission(j,obs(k,t+1,r),r);
          }
          if(tmp > neginf){
            sumtmp = logSumExp(sumtmp,tmp);
          }
        }
        beta(i,t,k) = sumtmp;        
      }
    }
  }
}


void internalBackward(const arma::mat& transition, const arma::mat& emission, 
                       const arma::Mat<int>& obs, arma::cube& beta) {  
  
  double sumtmp;
  double tmp;
  double neginf = -arma::math::inf();
  for(unsigned int k = 0; k < obs.n_rows; k++){
    beta.slice(k).col(beta.n_cols-1).fill(0.0);
    for(int t = (obs.n_cols-2); t >= 0; t--){  
      for(unsigned int i = 0; i < beta.n_rows; i++){
        sumtmp = neginf;
        for(unsigned int j = 0; j < beta.n_rows; j++){
          tmp = beta(j,t+1,k) + transition(i,j) + emission(j,obs(k,t+1));
          if(tmp > neginf){
            sumtmp = logSumExp(sumtmp,tmp);
          }
        }
        beta(i,t,k) = sumtmp;        
      }
    }
  }
}
