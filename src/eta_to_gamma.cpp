#include "eta_to_gamma.h"
#include "create_Q.h"
#include "sum_to_zero.h"

// [[Rcpp::export]]
arma::mat eta_to_gamma_mat(const arma::mat& eta) {
  arma::mat Q = create_Q(eta.n_rows + 1);
  return sum_to_zero(eta, Q);
}
// [[Rcpp::export]]
arma::cube eta_to_gamma_cube(const arma::cube& eta) {
  arma::mat Q = create_Q(eta.n_rows + 1);
  arma::cube gamma(eta.n_rows + 1, eta.n_cols, eta.n_slices);
  for (unsigned int i = 0; i < eta.n_slices; i++) {
    gamma.slice(i) = sum_to_zero(eta.slice(i), Q);
  }
  return gamma;
}
// [[Rcpp::export]]
arma::field<arma::mat> eta_to_gamma_mat_field(const arma::field<arma::mat>& eta) {
  unsigned int L = eta.n_elem;
  arma::field<arma::mat> gamma(L);
  for (unsigned int l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l));
  }
  return gamma;
}
// [[Rcpp::export]]
arma::field<arma::cube> eta_to_gamma_cube_field(const arma::field<arma::cube>& eta) {
  unsigned int L = eta.n_elem;
  arma::field<arma::cube> gamma(L);
  for (unsigned int l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l));
  }
  return gamma;
}

arma::mat eta_to_gamma(const arma::mat& eta) {
  arma::mat Q = create_Q(eta.n_rows + 1);
  return sum_to_zero(eta, Q);
}
arma::cube eta_to_gamma(const arma::cube& eta) {
  arma::mat Q = create_Q(eta.n_rows + 1);
  arma::cube gamma(eta.n_rows + 1, eta.n_cols, eta.n_slices);
  for (unsigned int i = 0; i < eta.n_slices; i++) {
    gamma.slice(i) = sum_to_zero(eta.slice(i), Q);
  }
  return gamma;
}
arma::field<arma::mat> eta_to_gamma(const arma::field<arma::mat>& eta) {
  unsigned int L = eta.n_elem;
  arma::field<arma::mat> gamma(L);
  for (unsigned int l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l));
  }
  return gamma;
}
arma::field<arma::cube> eta_to_gamma(const arma::field<arma::cube>& eta) {
  unsigned int L = eta.n_elem;
  arma::field<arma::cube> gamma(L);
  for (unsigned int l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l));
  }
  return gamma;
}

arma::mat eta_to_gamma(const arma::mat& eta, const arma::mat& Q) {
  return sum_to_zero(eta, Q);
}
arma::cube eta_to_gamma(const arma::cube& eta, const arma::mat& Q) {
  arma::cube gamma(eta.n_rows + 1, eta.n_cols, eta.n_slices);
  for (unsigned int i = 0; i < eta.n_slices; i++) {
    gamma.slice(i) = sum_to_zero(eta.slice(i), Q);
  }
  return gamma;
}
arma::field<arma::mat> eta_to_gamma(
    const arma::field<arma::mat>& eta, const arma::mat& Q) {
  unsigned int L = eta.n_elem;
  arma::field<arma::mat> gamma(L);
  for (unsigned int l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l), Q);
  }
  return gamma;
}
arma::field<arma::cube> eta_to_gamma(
    const arma::field<arma::cube>& eta, const arma::mat& Q) {
  unsigned int L = eta.n_elem;
  arma::field<arma::cube> gamma(L);
  for (unsigned int l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l), Q);
  }
  return gamma;
}