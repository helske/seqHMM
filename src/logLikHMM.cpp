
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logLikHMM(NumericVector transitionMatrix, NumericVector emissionArray, 
NumericVector initialProbs, IntegerVector obsArray) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p
  IntegerVector oDims = obsArray.attr("dim"); //k,n
  
  
  
  arma::colvec init(initialProbs.begin(),eDims[0]);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0]);
  arma::mat emission(emissionArray.begin(), eDims[0], eDims[1]);
  arma::Mat<int> obs(obsArray.begin(), oDims[0], oDims[1]);
  
  
  transition = log(transition); 
  emission = log(emission); 
  init = log(init); 
  
  double sumtmp;
  double tmp;
  double neginf = -arma::math::inf();
  
  arma::vec alpha(eDims[0]); //m,n,k
  arma::vec alphatmp(eDims[0]); //m,n,k
  
  double ll=0.0;
  
  for(int k = 0; k < oDims[0]; k++){
    
    
    for(int i=0; i < eDims[0]; i++){      
      alpha(i) = init(i)+emission(i,obs(k,0));
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
        alphatmp(i) = sumtmp + emission(i,obs(k,t));
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

