
#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List forwardbackward(NumericVector transitionMatrix, NumericVector emissionArray, 
  NumericVector initialProbs, IntegerVector obsArray, bool forwardonly, int threads = 1) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::colvec init(initialProbs.begin(),eDims[0],false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],true);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2],true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2],false);
 
  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k  
  arma::mat scales(oDims[1],oDims[0]); //n,k 
  
  internalForward(transition, emission, init, obs, alpha, scales, threads);
   
  if(forwardonly) {
    return List::create(Named("forward_probs") = wrap(alpha), Named("scaling_factors") = wrap(scales));
  } else {
    arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k
    internalBackward(transition, emission, obs, beta, scales, threads);
    return List::create(Named("forward_probs") = wrap(alpha),
      Named("backward_probs") = wrap(beta), Named("scaling_factors") = wrap(scales));
  }
  
}
