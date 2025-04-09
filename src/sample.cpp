#include "sample.h"
#include <RcppArmadilloExtensions/sample.h>

arma::uvec sample(const arma::uvec& x, const arma::vec& probs){
  return Rcpp::RcppArmadillo::sample(x, 1, false, probs);
}
