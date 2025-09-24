// softmax
#include "softmax.h"
#include "logsumexp.h"

// [[Rcpp::export]]
arma::vec log_softmax(const arma::vec& x) {
  return x - logSumExp(x);
}

// This version is recommended in
// Pierre Blanchard, Desmond J Higham, Nicholas J Higham (2021). 
// Accurately computing the log-sum-exp and softmax functions, 
// IMA Journal of Numerical Analysis, 41, 4, 2311–2330
// for better numerical accuracy than exp(log_softmax)
// [[Rcpp::export]]
arma::vec softmax(const arma::vec& x) {
  double x_max = arma::max(x);
  arma::vec result = arma::exp(x - x_max);
  return result / arma::accu(result);
}