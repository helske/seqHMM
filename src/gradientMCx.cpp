#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector gradientMCx(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
IntegerVector obsArray,NumericVector rowSumsA,NumericVector rowSumsB_, NumericVector sumInit,
IntegerVector transNZ, IntegerVector emissNZ, IntegerVector initNZ, NumericVector expPsi,
NumericMatrix coefs, NumericMatrix X_, IntegerVector numberOfStates) { 
  

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
  
 int q = coefs.nrow();
  arma::mat coef(coefs.begin(),q,numberOfStates.size());
  coef.col(0) = 0.0;
  arma::mat X(X_.begin(),oDims[0],q);
  
  arma::mat lweights = exp(X*coef).t();
  if(!lweights.is_finite()){
    return -std::numeric_limits<double>::max();
  }
  lweights.each_row() /= sum(lweights,0);
  lweights = log(lweights); 

  init = log(init); 
    
  arma::vec initk(eDims[0],oDims[0]);
  
  for(int k = 0; k < oDims[0]; k++){    
    initk.col(k) = init + reparma(lweights.col(k),numberOfStates);
  }
  
  
  arma::mat transitionLog = log(transition); 
  arma::cube emissionLog = log(emission);
  
   arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k 
  
  internalForwardMCx(transitionLog, emissionLog, initk, obs, alpha);
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
   arma::vec grad(expPsi.size()+q*(numberOfStates.size()-1),arma::fill::zeros);
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
              gradRow(j) += exp(initk(i,k)+tmp+beta(i,0,k)-ll(k));
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
        dpsi = (dpsi-init(ind).t())*expPsi(countgrad)/sumInit(0); //0 is wrong  
        grad(countgrad) = arma::as_scalar(dpsi*gradRow);
        countgrad ++;
      }
    }
  }
//  // beta
//  for(int k = 0; k < oDims[0]; k++){ 
//    for(unsigned int j = 1; j < eDims[0]; j++){
//      double tmp = 0.0;
//      double tmp2 = 0.0;
//      for(int i = 0; i< eDims[0]; i++){
//        if(i!=j){
//          tmp += arma::as_scalar(exp(dot(coef.col(i),X.row(k))));
//          for(int r=0; r < oDims[2]; r++){
//          tmp2 += arma::as_scalar(exp(emissionLog(i,obs(k,0,r),r)+beta(i,0,k)-ll(k)));
//          }
//        }
//      }
//      for(int r=0; r < oDims[2]; r++){
//       tmp *= exp(emissionLog(j,obs(k,0,r),r)+beta(j,0,k)-ll(k));
//      }
//      grad.subvec(expPsi.size()+q*(j-1),expPsi.size()+q*j-1) += X.row(k).t()*
//      arma::as_scalar(exp(dot(coef.col(j),X.row(k)))/pow(sumInit,2)*(tmp - tmp2));
//    }
//  }
  
  return wrap(grad);
}