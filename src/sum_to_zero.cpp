#include "sum_to_zero.h"

arma::mat sum_to_zero(const arma::mat& x, const arma::mat& Q) {
  arma::mat y(x.n_rows + 1, x.n_cols);
  for (arma::uword i = 0; i < x.n_cols; ++i) {
    y.col(i) = Q * x.col(i);
  }
  return y;
}
arma::vec sum_to_zero(const arma::vec& x, const arma::mat& Q) {
  return Q * x;
}
