#include "seqHMM.h"

// [[Rcpp::export]]

NumericMatrix varcoef(NumericMatrix coefs, NumericMatrix X_, IntegerVector numberOfStates) {

  int q = coefs.nrow();
  arma::mat coef(coefs.begin(), q, numberOfStates.size());
  coef.col(0).zeros();
  arma::mat X(X_.begin(), X_.nrow(), q);
  arma::mat weights = exp(X * coef).t();
  arma::rowvec sumweights = sum(weights, 0);

  weights.each_row() /= sumweights;
  return wrap(arma::inv_sympd(-hCoef(weights, X)));
}
