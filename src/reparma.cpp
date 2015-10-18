#include <algorithm>
#include "seqHMM.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
arma::vec reparma(arma::vec x, IntegerVector y) {
  int n = y.size();
  arma::vec res(sum(y));
  int ind = 0;
  for (int i = 0; i < n; ++i) {
    std::fill(res.begin() + ind, res.begin() + ind + y[i], x[i]);
    ind += y[i];
  }
  return res;
}
