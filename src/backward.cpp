
#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector backward(NumericVector transitionMatrix, NumericVector emissionArray, 
  IntegerVector obsArray) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],true);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2],true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2],false);
  
  
  transition = log(transition); 
  emission = log(emission);
  
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k
  internalBackward(transition, emission, obs, beta);
  
  return wrap(beta);
}
