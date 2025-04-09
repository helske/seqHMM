#include "config.h"
// [[Rcpp::export]]
arma::mat fast_quantiles(const arma::mat& X, const arma::vec& probs) {
  return arma::quantile(X, probs, 1);
}
