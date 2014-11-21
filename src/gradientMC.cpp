#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector gradientMC(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
IntegerVector obsArray,NumericVector rowSumsA,NumericVector rowSumsB_, double sumInit,
IntegerVector transNZ, IntegerVector emissNZ, IntegerVector initNZ, NumericVector expPsi) { 
  
  NumericVector grad(expPsi.size());
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(),eDims[0],false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],false);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1],eDims[2],false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1],oDims[2],false); 
  arma::imat ANZ(transNZ.begin(),eDims[0],eDims[0],false);
  arma::icube BNZ(emissNZ.begin(), eDims[0], eDims[1]-1,eDims[2],false);
  arma::ivec INZ(initNZ.begin(), eDims[0],false);
  arma::mat rowSumsB(rowSumsB_.begin(),eDims[0],eDims[2],false);
  
  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k 
  
  arma::vec initLog = log(init);
  arma::mat transitionLog = log(transition); 
  arma::cube emissionLog = log(emission);
  
  internalForwardMC(transitionLog, emissionLog, initLog, obs, alpha);
  internalBackwardMC(transitionLog, emissionLog, obs, beta);     
  
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
  
  int countgrad=0;
  // transitionMatrix
  for(int i = 0; i < eDims[0]; i++){   
    arma::uvec ind = arma::find(ANZ.row(i));
    
    if(ind.n_elem>1){ 
      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
      for(unsigned int j = 0; j < ind.n_elem; j++){
        for(int k = 0; k < oDims[0]; k++){ 
          for(int t = 0; t < (oDims[1]-1); t++){  
            tmp = 0.0;
            for(int r=0; r < oDims[2]; r++){
              tmp += emissionLog(ind(j),obs(k,t+1,r),r);
            }
            gradRow(j) += exp(alpha(i,t,k)+tmp+beta(ind(j),t+1,k)-ll(k));
          }
        }
      }      
      
      for(unsigned int j = 0; j < ind.n_elem; j++){
        if(ANZ(i,ind(j))!=2){
          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
          dpsi(j) = 1.0;    
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
  for(int r=0; r < oDims[2]; r++){
    for(int i = 0; i < eDims[0]; i++){   
      arma::uvec ind = arma::find(BNZ.slice(r).row(i)); 
      
      if(ind.n_elem>1){ 
        arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
        for(unsigned int j = 0; j < ind.n_elem; j++){
          for(int k = 0; k < oDims[0]; k++){               
            if(obs(k,0,r)==ind(j)){   
               tmp = 0.0;
              for(int r2=0; r2 < oDims[2]; r2++){
                if(r2!=r){
                  tmp += emissionLog(i,obs(k,0,r2),r2);
                }
              }
              gradRow(j) += exp(initLog(i)+tmp+beta(i,0,k)-ll(k));
            }       
            for(int t = 0; t < (oDims[1]-1); t++){              
              if(obs(k,t+1,r)==ind(j)){
                tmp = 0.0;
                for(int r2=0; r2 < oDims[2]; r2++){
                  if(r2!=r){
                    tmp += emissionLog(i,obs(k,t+1,r2),r2);
                  }
                }
                gradRow(j) += arma::accu(exp(alpha.slice(k).col(t)+tmp+transitionLog.col(i)+beta(i,t+1,k)-ll(k)));
              }          
            }
          }
        }
        for(unsigned int j = 0; j < ind.n_elem; j++){
          if(BNZ(i,ind(j),r)!=2){
            arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
            dpsi(j) = 1.0;    
            arma::uvec indi(1);
            indi(0)=i;
            dpsi = (dpsi-emission.slice(r).submat(indi,ind))*expPsi(countgrad)/rowSumsB(i,r);  
            grad(countgrad) = arma::as_scalar(dpsi*gradRow);
            countgrad ++;
          }
        }
      }
    }
    
  }
  // InitProbs  
  arma::uvec ind = arma::find(INZ);    
  if(ind.n_elem>1){ 
    arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
    for(unsigned int j = 0; j < ind.n_elem; j++){
      for(int k = 0; k < oDims[0]; k++){ 
        tmp = 0.0;
        for(int r=0; r < oDims[2]; r++){
          tmp += emissionLog(j,obs(k,0,r),r);
        }
        gradRow(j) += exp(tmp+beta(j,0,k)-ll(k));
      }
    }
    for(unsigned int j = 0; j < ind.n_elem; j++){
      if(INZ(ind(j))!=2){
        arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
        dpsi(j) = 1.0;    
        dpsi = (dpsi-init(ind).t())*expPsi(countgrad)/sumInit;  
        grad(countgrad) = arma::as_scalar(dpsi*gradRow);
        countgrad ++;
      }
    }
  }
  
  return grad;
}