
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector backward(NumericVector transitionMatrix, NumericVector emissionArray, 
		IntegerVector obsArray) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],true);
  arma::mat emission(emissionArray.begin(), eDims[0], eDims[1],true);
  arma::Mat<int> obs(obsArray.begin(), oDims[0], oDims[1],false);
  
  
  transition = log(transition); 
  emission = log(emission);
  
  double sumtmp;
  double tmp;
  double neginf = -arma::math::inf();
  
  for(int k = 0; k < oDims[0]; k++){    
    beta.slice(k).col(oDims[1]-1).fill(0.0);
    for(int t = oDims[1]-2; t >= 0; t--){  
      for(int i = 0; i < eDims[0]; i++){
        sumtmp = neginf;
        for(int j = 0; j < eDims[0]; j++){
          tmp = beta(j,t+1,k) + transition(i,j) + emission(j,obs(k,t+1));
          
         
          if(tmp > neginf){
            sumtmp = logSumExp(sumtmp,tmp);
          }
        }
        beta(i,t,k) = sumtmp;        
      }
    }
  }
  
  return wrap(beta);
}
