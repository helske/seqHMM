#include "sample.h"
#include <RcppArmadilloExtensions/sample.h>

arma::uword sample(const arma::uvec& x, const arma::rowvec& probs){
  return arma::as_scalar(Rcpp::RcppArmadillo::sample(x, 1, false, probs.t()));
}
