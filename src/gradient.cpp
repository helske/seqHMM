#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List gradient(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
IntegerVector obsArray,NumericVector rowSumsA,NumericVector rowSumsB, double sumInit,
IntegerVector transNZ, IntegerVector emissNZ, IntegerVector initNZ, NumericVector expPsi) { 
    
  arma::vec grad(expPsi.size(),arma::fill::zeros);
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(),eDims[0],false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],false);
  arma::mat emission(emissionArray.begin(), eDims[0], eDims[1],false);
  arma::imat obs(obsArray.begin(), oDims[0], oDims[1],false); 
  arma::imat ANZ(transNZ.begin(),eDims[0],eDims[0],false);
  arma::imat BNZ(emissNZ.begin(), eDims[0], eDims[1]-1,false);
  arma::ivec INZ(initNZ.begin(), eDims[0],false);
  
  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::mat scales(oDims[1],oDims[0]);
  
  internalForward(transition, emission, init, obs, alpha, scales);
  internalBackward(transition, emission, obs, beta, scales);     
  
 
  int countgrad=0;
  // transitionMatrix
  for(unsigned int i = 0; i < eDims[0]; i++){   
    arma::uvec ind = arma::find(ANZ.row(i));
    
    if(ind.n_elem>1){ 
      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
      for(unsigned int j = 0; j < ind.n_elem; j++){
        for(int k = 0; k < oDims[0]; k++){ 
          for(int t = 0; t < (oDims[1]-1); t++){  
            gradRow(j) += alpha(i,t,k)*emission(ind(j),obs(k,t+1))*beta(ind(j),t+1,k);
          }
        }
      }      
      
      for(unsigned int k = 0; k < ind.n_elem; k++){
        if(ANZ(i,ind(k))!=2){
          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
          dpsi(k) = 1.0;    
          arma::uvec indi(1);
          indi(0)=i;
          dpsi = (dpsi-transition(indi,ind))*expPsi(countgrad)/rowSumsA(i);  
          grad(countgrad) = arma::as_scalar(dpsi*gradRow);
          countgrad ++;
        }
      }
    }
  }
  // emissionMatrix
  for(unsigned int i = 0; i < eDims[0]; i++){   
    arma::uvec ind = arma::find(BNZ.row(i));    
    if(ind.n_elem>1){ 
      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
      for(unsigned int j = 0; j < ind.n_elem; j++){
        for(int k = 0; k < oDims[0]; k++){ 
          if(obs(k,0)==ind(j)){
            gradRow(j) += init(i)*beta(i,0,k);
            }
          for(int t = 0; t < (oDims[1]-1); t++){
            if(obs(k,t+1)==ind(j)){
            gradRow(j) += arma::accu(alpha.slice(k).col(t)%transition.col(i))*beta(i,t+1,k);
            }
          }
        }
      }
      for(unsigned int k = 0; k < ind.n_elem; k++){
        if(BNZ(i,ind(k))!=2){
          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
          dpsi(k) = 1.0;    
          arma::uvec indi(1);
          indi(0)=i;
          dpsi = (dpsi-emission(indi,ind))*expPsi(countgrad)/rowSumsB(i);  
          grad(countgrad) = arma::as_scalar(dpsi*gradRow);
          countgrad ++;
        }
      }
    }
  }
  
  return List::create(Named("logLik")=-arma::accu(log(scales)),Named("gradient")=grad,Named("scales")=scales,
  Named("beta")=beta);
}