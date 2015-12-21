//Variance-covariance matrix of beta coefficient based on inverse of -Hessian
#include "seqHMM.h"

// [[Rcpp::export]]

NumericMatrix varcoef(const arma::mat& coef, const arma::mat& X) {

  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= sum(weights, 0);
  // use inv instead of faster inv_sympd as the latter produces error on valgrind
  return wrap(arma::inv(-hCoef(weights, X)));
}
