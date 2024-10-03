#include "sum_to_zero.h"

arma::mat sum_to_zero(const arma::mat& x) {
  return arma::join_cols(x, -arma::sum(x));
}
arma::vec sum_to_zero(const arma::vec& x) {
  arma::vec z(x.n_elem + 1);
  z.rows(0, x.n_elem - 1) = x;
  z(x.n_elem) = - arma::accu(x);
  return z;
}
// [[Rcpp::export]]
arma::mat eta_to_gamma_mat(const arma::mat& x) {
  return sum_to_zero(x);
}
// [[Rcpp::export]]
arma::cube eta_to_gamma_cube(const arma::cube& x) {
  arma::cube z(x.n_rows + 1, x.n_cols, x.n_slices);
  for (unsigned int i = 0; i < x.n_slices; i++) {
    z.slice(i) = sum_to_zero(x.slice(i));
  }
  return z;
}
