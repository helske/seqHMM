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
  for (arma::uword i = 0; i < eta.n_slices; i++) {
    gamma.slice(i) = sum_to_zero(eta.slice(i), Q);
  }
  return gamma;
}
// [[Rcpp::export]]
arma::field<arma::mat> eta_to_gamma_mat_field(const arma::field<arma::mat>& eta) {
  arma::uword L = eta.n_elem;
  arma::field<arma::mat> gamma(L);
  for (arma::uword l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l));
  }
  return gamma;
}
// [[Rcpp::export]]
arma::field<arma::cube> eta_to_gamma_cube_field(const arma::field<arma::cube>& eta) {
  arma::uword L = eta.n_elem;
  arma::field<arma::cube> gamma(L);
  for (arma::uword l = 0; l < L; l++) {
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
  for (arma::uword i = 0; i < eta.n_slices; i++) {
    gamma.slice(i) = sum_to_zero(eta.slice(i), Q);
  }
  return gamma;
}
arma::field<arma::mat> eta_to_gamma(const arma::field<arma::mat>& eta) {
  arma::uword L = eta.n_elem;
  arma::field<arma::mat> gamma(L);
  for (arma::uword l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l));
  }
  return gamma;
}
arma::field<arma::cube> eta_to_gamma(const arma::field<arma::cube>& eta) {
  arma::uword L = eta.n_elem;
  arma::field<arma::cube> gamma(L);
  for (arma::uword l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l));
  }
  return gamma;
}

arma::mat eta_to_gamma(const arma::mat& eta, const arma::mat& Q) {
  return sum_to_zero(eta, Q);
}
arma::cube eta_to_gamma(const arma::cube& eta, const arma::mat& Q) {
  arma::cube gamma(eta.n_rows + 1, eta.n_cols, eta.n_slices);
  for (arma::uword i = 0; i < eta.n_slices; i++) {
    gamma.slice(i) = sum_to_zero(eta.slice(i), Q);
  }
  return gamma;
}
arma::field<arma::mat> eta_to_gamma(
    const arma::field<arma::mat>& eta, const arma::mat& Q) {
  arma::uword L = eta.n_elem;
  arma::field<arma::mat> gamma(L);
  for (arma::uword l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l), Q);
  }
  return gamma;
}
arma::field<arma::cube> eta_to_gamma(
    const arma::field<arma::cube>& eta, const arma::mat& Q) {
  arma::uword L = eta.n_elem;
  arma::field<arma::cube> gamma(L);
  for (arma::uword l = 0; l < L; l++) {
    gamma(l) = eta_to_gamma(eta(l), Q);
  }
  return gamma;
}


arma::cube rho_to_phi(const arma::cube& rho, const arma::mat& Q) {
  arma::cube phi(rho.n_rows + 1, rho.n_cols, rho.n_slices + 1);
  if (rho.n_cols > 0) {
    phi.slice(0).zeros(); // first category of y as reference
    for (arma::uword i = 1; i < (rho.n_slices + 1); i++) {
      phi.slice(i) = sum_to_zero(rho.slice(i - 1), Q);
    }
  }
  return phi;
}
arma::field<arma::cube> rho_to_phi(
    const arma::field<arma::cube>& rho, const arma::mat& Q) {
  arma::uword S = rho.n_elem;
  arma::field<arma::cube> phi(S);
  for (arma::uword s = 0; s < S; s++) {
    phi(s) = rho_to_phi(rho(s), Q);
  }
  return phi;
}
// [[Rcpp::export]]
arma::field<arma::cube> rho_to_phi_field(
    const arma::field<arma::cube>& rho) {
  arma::mat Q = create_Q(rho(0).n_rows + 1);
  arma::uword S = rho.n_elem;
  arma::field<arma::cube> phi(S);
  for (arma::uword s = 0; s < S; s++) {
    phi(s) = rho_to_phi(rho(s), Q);
  }
  return phi;
}