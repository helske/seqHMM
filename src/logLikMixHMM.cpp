
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logLikMixHMM(NumericVector transitionMatrix, NumericVector emissionArray, 
NumericVector initialProbs, IntegerVector obsArray, NumericMatrix coefs, 
NumericMatrix X_) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p
  IntegerVector oDims = obsArray.attr("dim"); //k,n  
  
  arma::colvec init(initialProbs.begin(),eDims[0]);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0]);
  arma::mat emission(emissionArray.begin(), eDims[0], eDims[1]);
  arma::Mat<int> obs(obsArray.begin(), oDims[0], oDims[1]);
  
  double tmp;
  arma::vec alpha(eDims[0]); //m,n,k
  arma::vec alphatmp(eDims[0]); //m,n,k
  
  double ll=0.0;
  
  int q = coefs.nrow();
  arma::mat coef(coefs.begin(),q,eDims[0]);
  coef.col(0) = 0.0;
  arma::mat X(X_.begin(),oDims[0],q);
  
  arma::mat lweights = exp(X*coef).t(); 
  lweights.each_row() /= sum(lweights,0);
  lweights = log(lweights); 
  
  transition = log(transition); 
  emission = log(emission); 
  init = log(init); 
  
  double sumtmp;  
  double neginf = -arma::math::inf();  
  
  for(int k = 0; k < oDims[0]; k++){    
    
    alpha = init+emission.col(obs(k,0));
    
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
    ll += tmp+lweights(k);
  }
  
  return ll;
}

