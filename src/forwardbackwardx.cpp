
#include "seqHMM.h"
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List forwardbackwardx(NumericVector transitionMatrix, NumericVector emissionArray, 
  NumericVector initialProbs, IntegerVector obsArray, NumericMatrix coefs, 
  NumericMatrix X_, IntegerVector numberOfStates, bool forwardonly) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::colvec init(initialProbs.begin(),eDims[0],true);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],true);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2],true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2],false);
  
  int q = coefs.nrow();
  arma::mat coef(coefs.begin(),q,coefs.ncol());
  coef.col(0).zeros();
  arma::mat X(X_.begin(),oDims[0],q);
  arma::mat lweights = exp(X*coef).t();
  lweights.each_row() /= arma::sum(lweights,0);
  
  arma::mat initk(eDims[0],oDims[0]);
  for(int k = 0; k < oDims[0]; k++){    
    initk.col(k) = init % reparma(lweights.col(k),numberOfStates);
  }
  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::mat scales(oDims[1],oDims[0]); //m,n,k
  internalForwardx(transition, emission, initk, obs, alpha, scales);
  
  if(forwardonly) {
    return List::create(Named("forward_probs") = wrap(alpha), Named("scaling_factors") = wrap(scales));
  } else {
    arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k
    internalBackward(transition, emission, obs, beta, scales);
    return List::create(Named("forward_probs") = wrap(alpha),
      Named("backward_probs") = wrap(beta), Named("scaling_factors") = wrap(scales));
  }
  
  return wrap(alpha);
}