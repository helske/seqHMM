//Variance-covariance matrix of beta coefficient based on inverse of -Hessian
#include "seqHMM.h"

// [[Rcpp::export]]

NumericMatrix varcoef(const arma::mat& coef, const arma::mat& X) {

  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= sum(weights, 0);
  return wrap(arma::inv_sympd(-hCoef(weights, X)));
}
