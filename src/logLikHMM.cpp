
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector logLikHMM(NumericVector transitionMatrix, NumericVector emissionArray, 
  NumericVector initialProbs, IntegerVector obsArray, int threads = 1) {  
  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  
  arma::colvec init(initialProbs.begin(),eDims[0], false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0], false);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  

  NumericVector ll(oDims[0]);  
  
#pragma omp parallel for num_threads(threads)
  for(int k = 0; k < oDims[0]; k++){    
    arma::vec alpha(eDims[0]);
    for(int i=0; i < eDims[0]; i++){      
      alpha(i) = init(i);
      for(int r = 0; r < oDims[2]; r++){
        alpha(i) *= emission(i,obs(k,0,r),r);
      }
    } 
    double tmp = sum(alpha);
    ll(k) = log(tmp);
    alpha /= tmp;
    
    arma::vec alphatmp(eDims[0]);
    
    for(int t = 1; t < oDims[1]; t++){  
      for(int i = 0; i < eDims[0]; i++){
        alphatmp(i) = arma::dot(transition.col(i), alpha);
        for(int r = 0; r < oDims[2]; r++){
          alphatmp(i) *= emission(i,obs(k,t,r),r);
        }
      }
      tmp = sum(alphatmp);
      ll(k) += log(tmp);
      alpha = alphatmp/tmp;
    }
  } 
  
  return ll;
}

