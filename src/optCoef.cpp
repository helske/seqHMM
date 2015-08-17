#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

void optCoef(arma::mat alpha, arma::mat beta, arma::vec& ll, arma::mat& coef, arma::mat& X, arma::mat& lweights, IntegerVector cumsumstate,
  IntegerVector numberOfStates) {
  
  arma::mat u = alpha + beta;
  for(int k = 0; k < alpha.n_cols; k++){
    u.col(k) -= ll(k);
  }
  u = exp(u);
  
  
  arma::mat usums(X.n_rows,numberOfStates.size());
  usums.zeros();
  
  for(int i = 0; i < numberOfStates.size(); i++){
    for(int k = 0; k < u.n_cols; k++){
      usums(k,i) += arma::accu(u.submat(cumsumstate(i)-numberOfStates(i),k,cumsumstate(i)-1,k));
    }
  }
  int p = X.n_cols;
  arma::vec tmp(p * (lweights.n_rows - 1));
  arma::mat coefnew(coef.n_rows,coef.n_cols - 1);
  int iter = 0;
  double change = 1.0;
  while((change>1e-10) & (iter<1000)){   
    lweights = exp(X*coef).t();
    lweights.each_row() /= sum(lweights,0);
    tmp = solve(hCoef(usums, lweights, X), gCoef(usums, lweights, X));
    for(int i = 0; i < (lweights.n_rows - 1); i++){
      coefnew.col(i) = coef.col(i + 1) - tmp.subvec(i * p, (i + 1) * p - 1);
    }
    change = arma::accu(arma::abs(coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) - coefnew))/coefnew.n_elem;
    coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) = coefnew;
    iter++;
    //Rcout<<iter<<" "<<change<<std::endl;
  }
  lweights = exp(X*coef).t();
  lweights.each_row() /= sum(lweights,0);
  lweights = log(lweights); 
}


arma::vec gCoef(arma::mat& usums, arma::mat& lweights, arma::mat& X) {
  
  int p = X.n_cols;
  arma::vec grad(p * (lweights.n_rows - 1)); 
  grad.zeros();
  for(int j = 0; j < (lweights.n_rows - 1); j++){
    for(int i = 0; i < X.n_rows; i++){
      grad.subvec(j*p,(j+1)*p-1) += X.row(i).t()*(usums(i,j+1) - lweights(j+1,i));
    }
  }
  
  return grad;
}

arma::mat hCoef(arma::mat& usums, arma::mat& lweights, arma::mat& X) {
  
  int p = X.n_cols;
  arma::mat hess(p * (lweights.n_rows - 1), p * (lweights.n_rows - 1));
  hess.zeros();
  for(int j = 0; j < (lweights.n_rows - 1); j++){
    for(int k = 0; k < (lweights.n_rows - 1); k++){
      for(int i = 0; i < X.n_rows; i++){ //X.n_rows
        if(j!=k){
          hess.submat(j*p,k*p,(j+1)*p-1,(k+1)*p-1) += X.row(i).t()*X.row(i)*lweights(j+1,i)*lweights(k+1,i);
        } else{
          hess.submat(j*p,j*p,(j+1)*p-1,(j+1)*p-1) -= X.row(i).t()*X.row(i)*lweights(j+1,i)*(1.0-lweights(j+1,i));
        
        }
        
      }
    }
  }
  return hess;
}


