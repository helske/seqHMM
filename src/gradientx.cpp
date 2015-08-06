#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector gradientx(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
IntegerVector obsArray,NumericVector rowSumsA,NumericVector rowSumsB, NumericVector sumInit,
IntegerVector transNZ, IntegerVector emissNZ, IntegerVector initNZ, NumericVector expPsi,
NumericMatrix coefs, NumericMatrix X_, IntegerVector numberOfStates) { 
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(),eDims[0],false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],false);
  arma::mat emission(emissionArray.begin(), eDims[0], eDims[1],false);
  arma::imat obs(obsArray.begin(), oDims[0], oDims[1],false); 
  arma::imat ANZ(transNZ.begin(),eDims[0],eDims[0],false);
  arma::imat BNZ(emissNZ.begin(), eDims[0], eDims[1]-1,false);
  arma::ivec INZ(initNZ.begin(), eDims[0],false);
  
  
  int q = coefs.nrow();
  arma::mat coef(coefs.begin(),q,numberOfStates.size());
  coef.col(0).zeros();
  arma::mat X(X_.begin(),oDims[0],q);
  
  arma::vec grad(expPsi.size()+q*(numberOfStates.size()-1),arma::fill::zeros);
  
  arma::mat lweights = exp(X*coef).t();
  if(!lweights.is_finite()){
    grad.fill(-std::numeric_limits<double>::max())
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
  arma::mat emissionLog = log(emission);
  
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
  
  
  
  int countgrad=0;
  // transitionMatrix
  for(int i = 0; i < eDims[0]; i++){   
    arma::uvec ind = arma::find(ANZ.row(i));
    
    if(ind.n_elem>1){ 
      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
      for(unsigned int j = 0; j < ind.n_elem; j++){
        for(int k = 0; k < oDims[0]; k++){ 
          for(int t = 0; t < (oDims[1]-1); t++){  
            gradRow(j) += exp(alpha(i,t,k)+emissionLog(ind(j),obs(k,t+1))+beta(ind(j),t+1,k)-ll(k));
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
  for(unsigned int i = 0; i < eDims[0]; i++){   
    arma::uvec ind = arma::find(BNZ.row(i));    
    if(ind.n_elem>1){ 
      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
      for(unsigned int j = 0; j < ind.n_elem; j++){
        for(int k = 0; k < oDims[0]; k++){ 
          if(obs(k,0)==ind(j)){
            gradRow(j) += exp(initk(i,k)+beta(i,0,k)-ll(k)); //exp(initLog(i)+beta(i,0,k)-ll(k)); 
          }
          for(int t = 0; t < (oDims[1]-1); t++){
            if(obs(k,t+1)==ind(j)){
              gradRow(j) += arma::accu(exp(alpha.slice(k).col(t)+transitionLog.col(i)+beta(i,t+1,k)-ll(k)));
            }
          }
        }
      }
      for(unsigned int j = 0; j < ind.n_elem; j++){
        if(BNZ(i,ind(j))!=2){
          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
          dpsi(j) = 1.0;    
          arma::uvec indi(1);
          indi(0)=i;
          dpsi = (dpsi-emission(indi,ind))*expPsi(countgrad)/rowSumsB(i);  
          grad(countgrad) = arma::as_scalar(dpsi*gradRow);
          countgrad ++;
        }
      }
    }
  }
  // InitProbs  
  IntegerVector cumsumstate = cumsum(numberOfStates);
    
  for(unsigned int jj = 0; jj < numberOfStates.size(); jj++){
    arma::uvec ind = arma::find(INZ.subvec(cumsumstate(jj)-numberOfStates(jj),
    cumsumstate(jj)-1));
    if(ind.n_elem>1){        
      arma::vec gradRow(numberOfStates(jj),arma::fill::zeros);  
      for(unsigned int j = 0; j < numberOfStates(jj); j++){  
        for(int k = 0; k < oDims[0]; k++){ 
          gradRow(j) += exp(emissionLog(cumsumstate(jj)-numberOfStates(jj)+j,obs(k,0))+
          beta(cumsumstate(jj)-numberOfStates(jj)+j,0,k)-ll(k)+lweights(jj,k));
        }
      }
      
      for(unsigned int j = 0; j < ind.n_elem; j++){
        if(INZ(cumsumstate(jj)-numberOfStates(jj)+ind(j))!=2){
          arma::rowvec dpsi(numberOfStates(jj),arma::fill::zeros);
          dpsi(ind(j)) = 1.0;    
          dpsi = (dpsi-init(arma::span(cumsumstate(jj)-numberOfStates(jj),cumsumstate(jj)-1)).t())*expPsi(countgrad)/sumInit(jj);
          grad(countgrad) = arma::as_scalar(dpsi*gradRow);          
          countgrad ++;
          
        }
      }
    }
  }
  // beta
  double tmp2;
  for(unsigned int jj = 1; jj < numberOfStates.size(); jj++){
    for(int k = 0; k < oDims[0]; k++){
      tmp2 = 0.0;
      for(unsigned int j = 0; j< numberOfStates.size(); j++){
        if(j!=jj){
          tmp2 += exp(dot(coef.col(j),X.row(k)));
        }
      }
      for(unsigned int j = 0; j < eDims[0]; j++){                       
        if(j>=(cumsumstate(jj)-numberOfStates(jj)) & (j<cumsumstate(jj))){
          grad.subvec(expPsi.size()+q*(jj-1),expPsi.size()+q*jj-1) += 
          exp(emissionLog(j,obs(k,0))+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*tmp2/sumweights(k);
        } else {
          grad.subvec(expPsi.size()+q*(jj-1),expPsi.size()+q*jj-1) -= 
          exp(emissionLog(j,obs(k,0))+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*
          exp(dot(coef.col(jj),X.row(k)))/sumweights(k);
        }
      }
    }
  }
  
  return wrap(grad);
}