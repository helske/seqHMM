
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logLikMCHMMx(NumericMatrix transitionMatrix, NumericVector emissionArray,
IntegerVector obsArray, NumericMatrix coef_, NumericMatrix X_) {  
  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0]); 
  transition = log(transition); 
 
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2]);
  arma::Cube<int> obs(obsArray.begin(), oDims[0], oDims[1], oDims[2]);
   emission = log(emission); 
   
  int q = coef_.nrow();
  arma::mat coef(coef_.begin(),q,eDims[0]);
  arma::mat X(X_.begin(),oDims[0],q);
  
  arma::mat init = exp(X*coef).t(); // m x k
  init.each_row() /= sum(init,0);
  init = log(init); 
  
  arma::vec alpha(eDims[0]); //m,n,k
  arma::vec alphatmp(eDims[0]); //m,n,k  
  double tmp;
  double ll=0.0;
  
  

  
  double sumtmp; 
  double neginf = -arma::math::inf();    
  for(int k = 0; k < oDims[0]; k++){    
    
    for(int i=0; i < eDims[0]; i++){      
      alpha(i) = init(i,k);
      for(int r = 0; r < oDims[2]; r++){
        alpha(i) += emission(i,obs(k,0,r),r);
      }
    }    
    
    for(int t = 1; t < oDims[1]; t++){  
      for(int i = 0; i < eDims[0]; i++){
        sumtmp = neginf;
        for(int j = 0; j < eDims[0]; j++){
          tmp = alpha(j) + transition(j,i);
          if(tmp > neginf){
            sumtmp = logSumExp(sumtmp,tmp);
          }
        }
        
        for(int r = 0; r < oDims[2]; r++){
          sumtmp += emission(i,obs(k,t,r),r);
        }
        alphatmp(i) = sumtmp;
      }
      alpha = alphatmp;
      
    }
    
    
    tmp = neginf;
    for(int i = 0; i < eDims[0]; i++){
      if(alpha(i)>neginf){
        tmp = logSumExp(alpha(i),tmp); 
      }
    }
    ll += tmp;    
    
  }
  
  
  
  return ll;
}

