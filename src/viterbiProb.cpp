
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double viterbiProb(NumericVector transitionMatrix, NumericVector emissionArray, 
NumericVector initialProbs, IntegerVector obsArray) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(), eDims[0], true);
  arma::mat transition(transitionMatrix.begin(), eDims[0], eDims[0], true);
  arma::mat emission(emissionArray.begin(), eDims[0], eDims[1], true);
  arma::Mat<int> obs(obsArray.begin(), oDims[0], oDims[1], true);
  
  init = log(init);
  transition = log(transition);
  emission = log(emission);
  
  double logp = 0.0;
  
  arma::vec deltaold(eDims[0]);
  arma::vec deltanew(eDims[0]);
  arma::uword phi;
  
  
  for(int k=0; k<oDims[0]; k++){
    
    deltaold = init + emission.col(obs(k,0));
    
    
    for(int t=1; t<oDims[1]; t++){
      for(int j=0; j<eDims[0]; j++){
        (deltaold+transition.col(j)).max(phi);
        deltanew(j) = deltaold(phi) + transition(phi,j) + emission(j,obs(k,t));
      }
      deltaold = deltanew;
    }
    
    logp += deltanew.max();
  }
  
  return logp;
}
