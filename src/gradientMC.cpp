#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector gradientMC(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
  IntegerVector obsArray,NumericVector rowSumsA,NumericVector rowSumsB_, double sumInit,
  IntegerVector transNZ, IntegerVector emissNZ, IntegerVector initNZ, NumericVector expPsi) { 
  
  
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
  arma::vec grad(expPsi.size());
  
  // transitionMatrix
  arma::vec gradArow(eDims[0]);
  arma::mat gradA(eDims[0],eDims[0]);
  for(int i = 0; i < eDims[0]; i++){
    arma::uvec ind = arma::find(ANZ.row(i));
    if(ind.n_elem>1){ 
      gradArow.zeros();
      gradA.zeros();
      gradA.diag() += arma::accu(transition.row(i));
      gradA.each_row() -= transition.row(i);
      gradA.each_row() /= arma::sum(transition,0);
      Rcout<<"i "<<i<<gradA<<std::endl;
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
      Rcout<<"gradArow "<<gradArow<<std::endl;
      gradArow = gradA * gradArow;
      grad.subvec(countgrad,countgrad+ind.n_elem-1) = gradArow.rows(ind);
      countgrad += ind.n_elem;
    }
  }
 
  // emissionMatrix
  arma::vec gradBrow(eDims[1]);
  arma::mat gradB(eDims[1],eDims[1]);
  for(int r=0; r < oDims[2]; r++){
    for(int i = 0; i < eDims[0]; i++){
      arma::uvec ind = arma::find(BNZ.slice(r).row(i));
      if(ind.n_elem>1){
        gradBrow.zeros();
        gradB.zeros();
        gradB.diag() += arma::accu(emission.slice(r).row(i));
        gradB.each_row() -= emission.slice(r).row(i);
        gradB.each_row() /= arma::sum(emission.slice(r).row(i),0);
        
        for(int k = 0; k < oDims[0]; k++){
          for(int j = 0; j < eDims[0]; j++){
            if(obs(k,0,r) == j){
              tmp = 0.0;
              for(int r2=0; r2 < oDims[2]; r2++){
                if(r2!=r){
                  tmp += emissionLog(i,obs(k,0,r2),r2);
                }
              }
              gradBrow(j) += exp(initLog(i) + tmp + beta(i,0,k) - ll(k));
            }
          }
          for(int t = 0; t < (oDims[1]-1); t++){ 
            for(int j = 0; j < eDims[0]; j++){
              if(obs(k,t+1,r) == j){
                tmp = 0.0;
                for(int r = 0; r < oDims[2]; r++){
                  tmp += emissionLog(j,obs(k,t+1,r),r);
                }
                gradBrow(j) += arma::accu(exp(alpha.slice(k).col(t) + tmp + transitionLog.col(i) + beta(i,t+1,k) - ll(k)));
              }
            }
          }
          gradBrow = gradB * gradBrow;
          
          grad.subvec(countgrad,countgrad+ind.n_elem-1) = gradBrow.rows(ind);
          countgrad += ind.n_elem;
        }
      }
    }
  }
  // InitProbs
  arma::uvec ind = arma::find(INZ); 
  if(ind.n_elem>1){
    arma::vec gradIrow(eDims[0]);
    arma::mat gradI(eDims[0],eDims[0]);
    
    
    gradIrow.zeros();
    gradI.zeros();
    gradI.diag() += arma::accu(init);
    gradI.each_row() -= init.t();
    tmp = as_scalar(arma::sum(init,0));
    for(int i = 0; i < eDims[0]; i++){
      gradI.row(i) /= tmp;
    }
    for(int k = 0; k < oDims[0]; k++){
      for(int j = 0; j < eDims[0]; j++){ 
        tmp = 0.0;
        for(int r = 0; r < oDims[2]; r++){
          tmp += emissionLog(j,obs(k,0,r),r);
        }
        gradIrow(j) += exp(tmp + beta(j,0,k) - ll(k)); 
      }
    }
    gradIrow = gradI * gradIrow;
    grad.subvec(countgrad,countgrad+ind.n_elem-1) = gradIrow.rows(ind);
    countgrad += ind.n_elem;
  }
  return wrap(grad);
}