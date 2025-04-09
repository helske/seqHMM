#include "create_Q.h"
// Creates matrix Q used for sum-to-zero transformation with orthogonal columns
// by transforming rbind(diag(S - 1), -1) using QR decomposition.
// Function to compute Givens rotation
// Simplified version of Algorithm 1.1 in 
// Hammarling, S., & Lucas, C. (2008). 
// Updating the QR factorization and the least squares problem.
// MIMS EPrint: 2008.111
// This takes into account the known structure of Q and R, i.e b != 0 always
arma::vec2 givens(const double a, const double b) {
  double c, s;
  if (std::abs(b) >= std::abs(a)) {
    double t = -a / b;
    s = 1 / std::sqrt(1 + t * t);
    c = s * t;
  } else {
    double t = -b / a;
    c = 1 / std::sqrt(1 + t * t);
    s = c * t;
  }
  return {c, s}; // Return c and s as a vector
}
// Simplified version of Algorithm 2.6 in 
// Hammarling, S., & Lucas, C. (2008).
// R and Q are identity matrices
arma::mat compute_cs(const arma::uword n) {
  arma::vec u(n, arma::fill::ones);
  u = -u;
  arma::mat R(n, n, arma::fill::eye);
  arma::mat cs(2, n);
  for (arma::uword j = 0; j < n; ++j) {
    cs.col(j) = givens(R(j, j), u(j));
    R(j, j) = cs(0, j) * R(j, j) - cs(1, j) * u(j);
    if (j < n - 1) {
      arma::vec t1 = R.col(j).rows(j + 1, n - 1);
      arma::vec t2 = u.rows(j + 1, n - 1);
      R.col(j).rows(j + 1, n - 1) = cs(0, j) * t1 - cs(1, j) * t2;
      u.rows(j + 1, n - 1) = cs(1, j) * t1 + cs(0, j) * t2;
    }
  }
  return cs;
}
arma::field<arma::mat> create_Q(const arma::uvec n) {
  arma::field<arma::mat> Q(n.n_elem);
  for (arma::uword i = 0; i < n.n_elem; ++i) {
    Q(i) = create_Q(n(i));
  }
  return Q;
}
// same as
// arma::mat Q(n, n, arma::fill::eye);
// Q.row(n - 1).fill(-1);
// and QR of Q
// [[Rcpp::export]]
arma::mat create_Q(const arma::uword n) {
  arma::mat cs = compute_cs(n - 1);
  arma::mat Q(n, n, arma::fill::eye);
  arma::vec t1(n);
  arma::vec t2(n);
  for (arma::uword j = 0; j < n - 1; ++j) {
    t1 = Q.col(j);
    t2 = Q.col(n - 1);
    Q.col(j) = cs(0, j) * t1 - cs(1, j) * t2;
    Q.col(n - 1) = cs(1, j) * t1 + cs(0, j) * t2;
  }
  return Q.cols(0, n - 2);
}
