//Variance-covariance matrix of gamma coefficients based on inverse of -Hessian
#include "optcoef.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix varcoef(const arma::mat& coef, const arma::mat& X) {

  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= sum(weights, 0);
  // use inv instead of faster inv_sympd as the latter produces error on valgrind
  return Rcpp::wrap(arma::inv(-hCoef(weights, X)));
}
