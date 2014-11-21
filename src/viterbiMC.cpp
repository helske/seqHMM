
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List viterbiMC(NumericVector transitionMatrix, NumericVector emissionArray, 
NumericVector initialProbs, IntegerVector obsArray) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(), eDims[0], true);
  arma::mat transition(transitionMatrix.begin(), eDims[0], eDims[0], true);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::Cube<int> obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], true);
  
  arma::umat q(oDims[0], oDims[1]);
  arma::vec logp(oDims[0]);
  
  arma::mat delta(eDims[0],oDims[1]);
  arma::umat phi(eDims[0],oDims[1]);
  
  
  for(int k=0; k<oDims[0]; k++){
    
    delta.col(0) = init;
    for(int r=0; r<eDims[2]; r++){
      delta.col(0) += emission.slice(r).col(obs(k,0,r));
    }   
    
    phi.col(0).zeros();
    
    for(int t=1; t<oDims[1]; t++){
      for(int j=0; j<eDims[0]; j++){
        (delta.col(t-1)+transition.col(j)).max(phi(j,t));
        delta(j,t) = delta(phi(j,t),t-1)+transition(phi(j,t),j);
        for(int r=0; r<eDims[2]; r++){
          delta(j,t) += emission(j,obs(k,t,r),r);
        }
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
