
#include "seqHMM.h"
using namespace Rcpp;

void internalForwardMC(const arma::mat& transition, const arma::cube& emission, 
const arma::vec& init, const arma::Cube<int>& obs, arma::cube& alpha, arma::mat& scales) {  
  
  double tmp;  
  arma::vec alphatmp(emission.n_rows);
  
  for(unsigned int k = 0; k < obs.n_rows; k++){    
    
    for(unsigned int i=0; i < emission.n_rows; i++){      
      alpha(i,0,k) = init(i);
      for(unsigned int r = 0; r < obs.n_slices; r++){
        alpha(i,0,k) *= emission(i,obs(k,0,r),r);
      }
    }
    scales(0,k) = 1.0/arma::accu(alpha.slice(k).col(0));
    alpha.slice(k).col(0) = alpha.slice(k).col(0)*scales(0,k);
    for(unsigned int t = 1; t < obs.n_cols; t++){  
      for(unsigned int i = 0; i < transition.n_rows; i++){    
        tmp = 0.0;
        for(unsigned int j = 0; j < transition.n_rows; j++){
          tmp += alpha(j,t-1,k)*transition(j,i);          
        }        
        for(unsigned int r = 0; r < obs.n_slices; r++){
          tmp *= emission(i,obs(k,t,r),r);
        }
        alphatmp(i) = tmp;
      }      
      scales(t,k) = 1.0/arma::accu(alphatmp);
      alpha.slice(k).col(t) = alphatmp*scales(t,k);
    }
    
  }
}


void internalForward(const arma::mat& transition, const arma::mat& emission, 
const arma::vec& init, const arma::Mat<int>& obs, arma::cube& alpha, arma::mat& scales) {  
  
  
  double tmp;
  arma::vec alphatmp(emission.n_rows);
  
  for(unsigned int k = 0; k < obs.n_rows; k++){    
    
    alpha.slice(k).col(0) = init%emission.col(obs(k,0)); 
    scales(0,k) = 1.0/arma::accu(alpha.slice(k).col(0));
    alpha.slice(k).col(0) = alpha.slice(k).col(0)*scales(0,k);
    for(unsigned int t = 1; t < obs.n_cols; t++){  
      for(unsigned int i = 0; i < transition.n_rows; i++){    
        tmp = 0.0;
        for(unsigned int j = 0; j < transition.n_rows; j++){
          tmp += alpha(j,t-1,k)*transition(j,i);          
        }                
        tmp *= emission(i,obs(k,t));
        
        alphatmp(i) = tmp;
      }      
      scales(t,k) = 1.0/arma::accu(alphatmp);
      alpha.slice(k).col(t) = alphatmp*scales(t,k);
    }    
    
  }
}

void internalForwardMCLog(const arma::mat& transition, const arma::cube& emission, 
const arma::vec& init, const arma::Cube<int>& obs, arma::cube& alpha) {  
  
  double sumtmp;
  double tmp;
  double neginf = -arma::math::inf();
  
  for(unsigned int k = 0; k < obs.n_rows; k++){      
    for(unsigned int i=0; i < emission.n_rows; i++){      
      alpha(i,0,k) = init(i);
      for(unsigned int r = 0; r < obs.n_slices; r++){
        alpha(i,0,k) += emission(i,obs(k,0,r),r);
      }
    }    
    for(unsigned int t = 1; t < obs.n_cols; t++){  
      for(unsigned int i = 0; i < transition.n_rows; i++){
        sumtmp = neginf;
        for(unsigned int j = 0; j < transition.n_rows; j++){
          tmp = alpha(j,t-1,k) + transition(j,i);
          if(tmp > neginf){
            sumtmp = logSumExp(sumtmp,tmp);
          }
        }
        alpha(i,t,k) = sumtmp;
        for(unsigned int r = 0; r < obs.n_slices; r++){
          alpha(i,t,k) += emission(i,obs(k,t,r),r);
        }
      }
    }    
  }
}


void internalForwardLog(const arma::mat& transition, const arma::mat& emission, 
const arma::vec& init, const arma::Mat<int>& obs, arma::cube& alpha) {  
  
  
  double sumtmp;
  double tmp;
  double neginf = -arma::math::inf();
  
  for(unsigned int k = 0; k < obs.n_rows; k++){
    
    for(unsigned int i=0; i < emission.n_rows; i++){
      alpha(i,0,k) = init(i) + emission(i,obs(k,0));
    }
    
    for(unsigned int t = 1; t < obs.n_cols; t++){  
      for(unsigned int i = 0; i < transition.n_rows; i++){
        sumtmp = neginf;
        for(unsigned int j = 0; j < transition.n_rows; j++){
          tmp = alpha(j,t-1,k) + transition(j,i);
          if(tmp > neginf){
            sumtmp = logSumExp(sumtmp,tmp);
          }
        }
        alpha(i,t,k) = sumtmp + emission(i,obs(k,t));          
      }
    }
  }
}