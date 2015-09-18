#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat optCoef(const arma::icube& obs, const arma::cube& emission, const arma::mat& initk, 
  const arma::cube& beta, const arma::vec& ll, arma::mat& coef, const arma::mat& X, 
  const IntegerVector cumsumstate, const IntegerVector numberOfStates, int trace) {
  
  arma::mat weights = exp(X*coef).t();
  weights.each_row() /= sum(weights,0);
  
  int p = X.n_cols;
  arma::vec tmpvec(p * (weights.n_rows - 1));
  arma::mat coefnew(coef.n_rows,coef.n_cols - 1);
  int iter = 0;
  double change = 1.0;
  while((change>1e-8) & (iter<100)){
    tmpvec = arma::solve(hCoef(weights, X), gCoef(obs, beta, emission, initk, weights, ll, X, cumsumstate, numberOfStates));
    for(int i = 0; i < (weights.n_rows - 1); i++){
      coefnew.col(i) = coef.col(i + 1) - tmpvec.subvec(i * p, (i + 1) * p - 1);
    }
    change = arma::accu(arma::abs(coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) - coefnew))/coefnew.n_elem;
    coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) = coefnew;
    iter++;
    if(trace==3){
      Rcout<<"coefficient optimization iter: "<< iter;
      Rcout<<" new coefficients: "<< std::endl<<coefnew<<std::endl;
      Rcout<<" relative change: "<<change<<std::endl;
    }
    weights = exp(X*coef).t();
    weights.each_row() /= sum(weights,0);
  }
  
  return(log(weights));
}


arma::vec gCoef(const arma::icube& obs, const arma::cube& beta, const arma::cube& emission, const arma::mat& initk,
  const arma::mat& weights, const arma::vec& ll, const arma::mat& X, const IntegerVector cumsumstate, const IntegerVector numberOfStates) {
  
  int q = X.n_cols;
  arma::vec grad(q * (weights.n_rows - 1) ,arma::fill::zeros);
  double tmp;
  for(unsigned int jj = 1; jj < numberOfStates.size(); jj++){
    for(int k = 0; k < obs.n_rows; k++){
      for(unsigned int j = 0; j < emission.n_rows; j++){                
        tmp = 0.0;
        for(int r=0; r < obs.n_slices; r++){
          tmp += emission(j,obs(k,0,r),r);
        }        
        if(j>=(cumsumstate(jj)-numberOfStates(jj)) & j<cumsumstate(jj)){
          grad.subvec(q*(jj-1),q*jj-1) += 
            exp(tmp+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*(1.0 - weights(jj,k)); 
        } else {
          grad.subvec(q*(jj-1),q*jj-1) -= 
            exp(tmp+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*weights(jj,k);
        }
      }
    }
  }
  
  return grad;
}

arma::mat hCoef(const arma::mat& weights, const arma::mat& X) {
  
  int p = X.n_cols;
  arma::mat hess(p * (weights.n_rows - 1), p * (weights.n_rows - 1));
  hess.zeros();
  for(int j = 0; j < (weights.n_rows - 1); j++){
    for(int k = 0; k < (weights.n_rows - 1); k++){
      for(int i = 0; i < X.n_rows; i++){
        if(j!=k){
          hess.submat(j*p,k*p,(j+1)*p-1,(k+1)*p-1) += X.row(i).t()*X.row(i)*weights(j+1,i)*weights(k+1,i);
        } else{
          hess.submat(j*p,j*p,(j+1)*p-1,(j+1)*p-1) -= X.row(i).t()*X.row(i)*weights(j+1,i)*(1.0-weights(j+1,i));
          
        }
        
      }
    }
  }
  return hess;
}


