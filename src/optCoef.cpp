#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

void optCoef(arma::mat alpha, arma::mat beta, arma::vec ll, arma::mat coef, arma::mat X, arma::mat lweights,IntegerVector cumsumstate,
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
 
     lweights = exp(X*coef).t();
     lweights.each_row() /= sum(lweights,0);
  //   int iter = 0;
  //   double change = 1.0;
  //   while((change>1e-5) & (iter<100)){   
  //     coef
  //   }
  Rcout<<hCoef(usums, lweights, X)<<std::endl;
  Rcout<<"4"<<std::endl;
  
  //   lweights = exp(X*coef).t();
  //   lweights.each_row() /= sum(lweights,0);
  //   
  //   lweights = log(lweights); 
}


arma::mat gCoef(arma::mat usums, arma::mat lweights, arma::mat X) {
  
  arma::mat grad(X.n_cols,lweights.n_rows-1); //6x2
  grad.zeros();
  for(int j=0; j< lweights.n_rows-1; j++){
    for(int i=0; i < X.n_rows; i++){
      grad.col(j) += X.row(i).t()*(usums(i,j+1) - lweights(j+1,i));
    }
  }
  
  return -grad;
}

arma::mat hCoef(arma::mat usums, arma::mat lweights, arma::mat X) {
  
  int p = X.n_cols;
  arma::mat hess(p*(lweights.n_rows-1),p*(lweights.n_rows-1));
  
  for(int j = 0; j < lweights.n_rows-1; j++){
    for(int k = 0; k < lweights.n_rows-1; k++){
      for(int i = 0; i < X.n_rows; i++){
        if(j!=k){
          hess.submat(j*p,k*p,(j+1)*p-1,(k+1)*p-1) += X.row(i).t()*X.row(i)*lweights(j+1,i)*lweights(k+1,i);
        } else{
          hess.submat(j*p,k*p,(j+1)*p-1,(k+1)*p-1) += X.row(i).t()*X.row(i)*lweights(j+1,i)*(1.0-lweights(j+1,i));
        }
        
      }
    }
  }
  return -hess;
}


