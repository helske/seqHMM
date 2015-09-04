#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector gradbeta(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
  IntegerVector obsArray, IntegerVector nSymbols,
  NumericMatrix coefs, NumericMatrix X_, IntegerVector numberOfStates) { 
  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(),eDims[0],false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],false);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1],eDims[2],false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1],oDims[2],false); 
  
  
  
  int q = coefs.nrow();
  arma::vec grad((numberOfStates.size()-1)*q ,arma::fill::zeros);
  arma::mat coef(coefs.begin(),q,numberOfStates.size());
  coef.col(0).zeros();
  arma::mat X(X_.begin(),oDims[0],q);
  arma::mat lweights = exp(X*coef).t();
  if(!lweights.is_finite()){
    grad.fill(-std::numeric_limits<double>::max());
    return wrap(grad);
  }
  arma::rowvec sumweights = sum(lweights,0);
  
  lweights.each_row() /= sumweights;
  lweights = log(lweights); 
  
  arma::vec initLog = log(init); 
  arma::mat initk(eDims[0],oDims[0]);
  for(int k = 0; k < oDims[0]; k++){    
    initk.col(k) = initLog + reparma(lweights.col(k),numberOfStates);
  }
  arma::mat transitionLog = log(transition); 
  arma::cube emissionLog = log(emission);
  
  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k 
  
  internalForwardx(transitionLog, emissionLog, initk, obs, alpha);
  internalBackward(transitionLog, emissionLog, obs, beta);     
  
  arma::vec ll(oDims[0]);
  
  double tmp=0.0;
  double neginf = -arma::math::inf();
  
  for(int k=0;k<oDims[0];k++){
    tmp =neginf;
    for(int i = 0; i < eDims[0]; i++){
      if(alpha(i,oDims[1]-1,k)>neginf){
        tmp = logSumExp(alpha(i,oDims[1]-1,k),tmp); 
      }
    }
    ll(k) = tmp;
  }
  
  IntegerVector cumsumstate = cumsum(numberOfStates);
  
  double tmp2;
  for(unsigned int jj = 1; jj < numberOfStates.size(); jj++){
    for(int k = 0; k < oDims[0]; k++){
      for(unsigned int j = 0; j < eDims[0]; j++){                
        tmp = 0.0;
        for(int r=0; r < oDims[2]; r++){
          tmp += emissionLog(j,obs(k,0,r),r);
        }        
        if(j>=(cumsumstate(jj)-numberOfStates(jj)) & j<cumsumstate(jj)){
          grad.subvec(q*(jj-1),q*jj-1) += 
            exp(tmp+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*(1.0 - exp(lweights(jj,k))); 
        } else {
          grad.subvec(q*(jj-1),q*jj-1) -= 
            exp(tmp+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*exp(lweights(jj,k));
        }
      }
    }
  }
  return wrap(-grad);
}