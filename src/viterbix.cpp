
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List viterbix(NumericVector transitionMatrix, NumericVector emissionArray, 
NumericVector initialProbs, IntegerVector obsArray, NumericMatrix coefs, 
NumericMatrix X_, IntegerVector numberOfStates) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p
  IntegerVector oDims = obsArray.attr("dim"); //k,n 
  
  arma::colvec init(initialProbs.begin(),eDims[0],false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],false);
  arma::mat emission(emissionArray.begin(), eDims[0], eDims[1],false);
  arma::imat obs(obsArray.begin(), oDims[0], oDims[1],false);
  arma::umat q(oDims[0], oDims[1]);
  arma::vec logp(oDims[0]);
  
  arma::mat delta(eDims[0],oDims[1]);
  arma::umat phi(eDims[0],oDims[1]);
  
  int qn = coefs.nrow();
  arma::mat coef(coefs.begin(),qn,numberOfStates.size());
 coef.col(0).zeros();
  arma::mat X(X_.begin(),oDims[0],qn);
  
  arma::mat lweights = exp(X*coef).t();
  lweights.each_row() /= sum(lweights,0);
  lweights = log(lweights);
  
  for(int k = 0; k < oDims[0]; k++){
    
    delta.col(0) = init + reparma(lweights.col(k),numberOfStates)+emission.col(obs(k,0));
    
    phi.col(0).zeros();
    
    for(int t=1; t<oDims[1]; t++){
      for(int j=0; j<eDims[0]; j++){
        (delta.col(t-1)+transition.col(j)).max(phi(j,t));
        delta(j,t) = delta(phi(j,t),t-1)+transition(phi(j,t),j)+emission(j,obs(k,t));
      }        
    }
    delta.col(oDims[1]-1).max(q(k,oDims[1]-1));
    for(int t=(oDims[1]-2); t>=0; t--){
      q(k,t) = phi(q(k,t+1),t+1);
    }
    logp(k) = delta.col(oDims[1]-1).max();
  }
  
  return List::create(Named("q") = wrap(q),Named("logp") = wrap(logp));
}