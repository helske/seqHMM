#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List objective(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
               IntegerVector obsArray, IntegerVector transNZ, IntegerVector emissNZ, 
               IntegerVector initNZ, IntegerVector nSymbols, int threads) { 
  
  
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
  
  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::mat scales(oDims[1],oDims[0]); //m,n,k
  
  internalForward(transition, emission, init, obs, alpha, scales, threads);
  if(!alpha.is_finite()){
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(),
                        Named("gradient") = wrap(grad));
  }
  internalBackward(transition, emission, obs, beta, scales, threads);     
  if(!beta.is_finite()){
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), 
                        Named("gradient") = wrap(grad));
  }
  arma::rowvec ll = arma::sum(log(scales));
  
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
            double tmp = 1.0;
            for(int r = 0; r < oDims[2]; r++){
              tmp *= emission(j,obs(k,t+1,r),r);
            }
            gradArow(j) += alpha(i,t,k) * tmp * beta(j,t+1,k) / scales(t+1,k); 
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
              double tmp = 1.0;
              for(int r2 = 0; r2 < oDims[2]; r2++){
                if(r2 != r){
                  tmp *= emission(i,obs(k,0,r2),r2);
                }
              }
              gradBrow(j) += init(i) * tmp * beta(i,0,k) / scales(0,k);
            }
            for(int t = 0; t < (oDims[1]-1); t++){ 
              if(obs(k,t+1,r) == j){
                double tmp = 1.0;
                for(int r2 = 0; r2 < oDims[2]; r2++){
                  if(r2 != r){
                    tmp *= emission(i,obs(k,t+1,r2),r2);
                  }
                }
                gradBrow(j) += arma::dot(alpha.slice(k).col(t),transition.col(i)) * tmp * beta(i,t+1,k) / scales(t+1,k);
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
        double tmp = 1.0;
        for(int r = 0; r < oDims[2]; r++){
          tmp *= emission(j,obs(k,0,r),r);
        }
        gradIrow(j) += tmp * beta(j,0,k) / scales(0,k); 
      }
    }
    gradIrow = gradI * gradIrow;
    grad.subvec(countgrad,countgrad+ind.n_elem-1) = gradIrow.rows(ind);
    countgrad += ind.n_elem;
  }
  return List::create(Named("objective") = -sum(ll), Named("gradient") = wrap(-grad));
}