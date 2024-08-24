// softmax
#include "logsumexp.h"
#include "softmax.h"

// [[Rcpp::export]]
arma::vec softmax(const arma::vec& x, const int logspace) {
  arma::vec result;
  if (logspace == 0) {
    double x_max = arma::max(x);
    result = arma::exp(x - x_max);
    result = result / sum(result);
  } else {
    double x_max = arma::max(x);
    result = x - x_max;
    result = result - logSumExp(result);
  }
  return result;
}

