
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

void viterbiForEM(arma::mat& transition, arma::mat& emission, arma::vec& init, arma::imat& obs, 
  arma::vec& logp, arma::umat& q) {  
  
  
  arma::mat delta(emission.n_rows,obs.n_cols);
  arma::umat phi(emission.n_rows,obs.n_cols);
  
  for(int k=0; k<obs.n_rows; k++){
    
    delta.col(0) = init+emission.col(obs(k,0));
    
    phi.col(0).zeros();
    
    for(int t=1; t<obs.n_cols; t++){
      for(int j=0; j<emission.n_rows; j++){
        (delta.col(t-1)+transition.col(j)).max(phi(j,t));
        delta(j,t) = delta(phi(j,t),t-1)+transition(phi(j,t),j)+emission(j,obs(k,t));
      }        
    }
    delta.col(obs.n_cols-1).max(q(k,obs.n_cols-1));
    for(int t=(obs.n_cols-2); t>=0; t--){
      q(k,t) = phi(q(k,t+1),t+1);
    }
    logp(k) = delta.col(obs.n_cols-1).max();
  }
  
}
