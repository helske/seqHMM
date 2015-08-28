#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector objectivex(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
  IntegerVector obsArray, IntegerVector transNZ, IntegerVector emissNZ, IntegerVector initNZ, IntegerVector nSymbols,
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
  
  arma::vec grad(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ),arma::fill::zeros);
  
  int q = coefs.nrow();
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
  
  int countgrad = 0;

  // transitionMatrix
  arma::vec gradArow(eDims[0]);
  arma::mat gradA(eDims[0],eDims[0]);
  for(int i = 0; i < eDims[0]; i++){
    arma::uvec ind = arma::find(ANZ.row(i));
    
    if(ind.n_elem>0){ 
      gradArow.zeros();
      gradA.eye();
      gradA.each_row() -= transition.row(i);
      gradA.each_col() %= transition.row(i).t();
      
      for(int k = 0; k < oDims[0]; k++){
        for(int t = 0; t < (oDims[1]-1); t++){
          for(int j = 0; j < eDims[0]; j++){ 
            tmp = 0.0;
            for(int r = 0; r < oDims[2]; r++){
              tmp += emissionLog(j,obs(k,t+1,r),r);
            }
            gradArow(j) += exp(alpha(i,t,k) + tmp + beta(j,t+1,k) - ll(k)); 
          }
          
        }
      }
      gradArow = gradA * gradArow;
      grad.subvec(countgrad,countgrad+ind.n_elem-1) = gradArow.rows(ind);
      countgrad += ind.n_elem;
    }
  }
  // emissionMatrix
  for(int r=0; r < oDims[2]; r++){
    arma::vec gradBrow(nSymbols[r]);
    arma::mat gradB(nSymbols[r],nSymbols[r]);
    for(int i = 0; i < eDims[0]; i++){
      arma::uvec ind = arma::find(BNZ.slice(r).row(i));
      if(ind.n_elem>0){
        gradBrow.zeros();
        gradB.eye();
        gradB.each_row() -= emission.slice(r).row(i).subvec(0,nSymbols[r]-1);
        gradB.each_col() %= emission.slice(r).row(i).subvec(0,nSymbols[r]-1).t();
        for(unsigned int j = 0; j < nSymbols[r]; j++){
          for(int k = 0; k < oDims[0]; k++){
            if(obs(k,0,r) == j){
              tmp = 0.0;
              for(int r2 = 0; r2 < oDims[2]; r2++){
                if(r2 != r){
                  tmp += emissionLog(i,obs(k,0,r2),r2);
                }
              }
              gradBrow(j) += exp(initk(i,k) + tmp + beta(i,0,k) - ll(k));
            }
            for(int t = 0; t < (oDims[1]-1); t++){ 
              if(obs(k,t+1,r) == j){
                tmp = 0.0;
                for(int r2 = 0; r2 < oDims[2]; r2++){
                  if(r2 != r){
                    tmp += emissionLog(i,obs(k,t+1,r2),r2);
                  }
                }
                gradBrow(j) += arma::accu(exp(alpha.slice(k).col(t) + tmp + transitionLog.col(i) + beta(i,t+1,k) - ll(k)));
              }
            }
          }
        }
        gradBrow = gradB * gradBrow;
        grad.subvec(countgrad,countgrad+ind.n_elem-1) = gradBrow.rows(ind);
        countgrad += ind.n_elem;
        
      }
    }
  }
  // InitProbs
  arma::uvec ind = arma::find(INZ); 
  if(ind.n_elem>0){
    arma::vec gradIrow(eDims[0]);
    arma::mat gradI(eDims[0],eDims[0]);
    
    gradIrow.zeros();
    gradI.zeros();
    gradI.eye();
    gradI.each_row() -= init.t();
    gradI.each_col() %= init;
    for(int k = 0; k < oDims[0]; k++){
      for(int j = 0; j < eDims[0]; j++){ 
        tmp = 0.0;
        for(int r = 0; r < oDims[2]; r++){
          tmp += emissionLog(j,obs(k,0,r),r);
        }
        gradIrow(j) += exp(tmp + beta(j,0,k) - ll(k) + lweights(j,k)); 
      }
    }
    gradIrow = gradI * gradIrow;
    grad.subvec(countgrad,countgrad+ind.n_elem-1) = gradIrow.rows(ind);
    countgrad += ind.n_elem;
  }
  
  
  // beta
  arma::vec gradBetarow(numberOfStates.size());
  arma::mat gradBeta(numberOfStates.size(),numberOfStates.size());
  IntegerVector cumsumstate = cumsum(numberOfStates);
  for(int k = 0; k < oDims[0]; k++){
    gradBetarow.zeros();
    gradBeta.zeros();
    gradBeta.each_col() -= lweights.col(k);
    gradBeta.each_row() %= lweights.col(k).t();
    gradBeta.diag() += 1.0;
    
    for(unsigned int i = 0; i < numberOfStates.size(); i++){
      gradBetarow(i) = arma::accu(exp(alpha.slice(k).col(oDims[1]-1).subvec(cumsumstate(i)-numberOfStates(i),cumsumstate(i)-1)));
    }
    
    gradBetarow = gradBeta * gradBetarow;
    grad.subvec(countgrad,grad.n_elem-1) += kron(gradBetarow, X.col(k));
  } 
  
  return wrap(grad);
}