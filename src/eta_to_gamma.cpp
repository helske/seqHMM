#include "eta_to_gamma.h"
#include "create_Q.h"
#include "sum_to_zero.h"

arma::mat eta_to_gamma(const arma::mat& eta, const arma::mat& Q) {
  return sum_to_zero(eta, Q);
}
arma::cube eta_to_gamma(const arma::cube& eta, const arma::mat& Q) {
  arma::cube gamma(eta.n_rows + 1, eta.n_cols, eta.n_slices);
  for (arma::uword i = 0; i < eta.n_slices; ++i) {
    gamma.slice(i) = sum_to_zero(eta.slice(i), Q);
  }
  return gamma;
}
arma::field<arma::mat> eta_to_gamma(const arma::field<arma::mat>& eta, const arma::mat& Q) {
  arma::uword D = eta.n_rows;
  arma::field<arma::mat> gamma(D);
  for (arma::uword d = 0; d < D; ++d) {
    gamma(d) = eta_to_gamma(eta(d), Q);
  }
  return gamma;
}
arma::field<arma::cube> eta_to_gamma(const arma::field<arma::cube>& eta, const arma::mat& Q) {
  arma::uword D = eta.n_rows;
  arma::field<arma::cube> gamma(D);
  for (arma::uword d = 0; d < D; ++d) {
    gamma(d) = eta_to_gamma(eta(d), Q);
  }
  return gamma;
}
// for eta_B
arma::field<arma::cube> eta_to_gamma(const arma::field<arma::cube>& eta, 
                                     const arma::field<arma::mat>& Q, 
                                     const arma::uword D) {
  if (D == 1) {
    arma::uword C = eta.n_rows;
    arma::field<arma::cube> gamma(C);
    for (arma::uword c = 0; c < C; ++c) {
      gamma(c) = eta_to_gamma(eta(c), Q(c));
    }
    return gamma;
  } else {
    arma::uword C = eta.n_cols;
    arma::field<arma::cube> gamma(D, C);
    for (arma::uword d = 0; d < D; ++d) {
      for (arma::uword c = 0; c < C; ++c) {
        gamma(d, c) = eta_to_gamma(eta(d, c), Q(c));
      }
    }
    return gamma;
  }
}

// for R
// for gamma_pi and gamma_omega
// [[Rcpp::export]]
arma::mat eta_to_gamma_mat(const arma::mat& eta) {
  arma::mat Q = create_Q(eta.n_rows + 1);
  return sum_to_zero(eta, Q);
}
// for gamma_pi in the mixture case
// [[Rcpp::export]]
arma::field<arma::mat> eta_to_gamma_mat_field(const arma::field<arma::mat>& eta) {
  arma::uword D = eta.n_elem;
  arma::field<arma::mat> gamma(D);
  arma::mat Q = create_Q(eta(0).n_rows + 1);
  for (arma::uword d = 0; d < D; ++d) {
    gamma(d) = eta_to_gamma(eta(d), Q);
  }
  return gamma;
}
// for gamma_A
// [[Rcpp::export]]
arma::cube eta_to_gamma_cube(const arma::cube& eta) {
  arma::mat Q = create_Q(eta.n_rows + 1);
  arma::cube gamma(eta.n_rows + 1, eta.n_cols, eta.n_slices);
  for (arma::uword i = 0; i < eta.n_slices; ++i) {
    gamma.slice(i) = eta_to_gamma(eta.slice(i), Q);
  }
  return gamma;
}
// for gamma_A in the mixture case and gamma_B in the non-mixture case
// [[Rcpp::export]]
arma::field<arma::cube> eta_to_gamma_cube_field(const arma::field<arma::cube>& eta) {
  arma::uword L = eta.n_elem;
  arma::field<arma::cube> gamma(L);
  for (arma::uword l = 0; l < L; ++l) {
    arma::mat Q = create_Q(eta(l).n_rows + 1);
    gamma(l) = eta_to_gamma(eta(l), Q);
  }
  return gamma;
}
// for gamma_B in the mixture case
// [[Rcpp::export]]
arma::field<arma::cube> eta_to_gamma_cube_2d_field(const Rcpp::List& eta) {
  arma::uword D = eta.size();
  Rcpp::List first_list = eta[0];
  arma::uword C = first_list.size();
  arma::field<arma::cube> gamma(D, C);
  for (arma::uword d = 0; d < D; ++d) {
    Rcpp::List eta_d = eta[d];
    for (arma::uword c = 0; c < C; ++c) {
      arma::cube eta_c = Rcpp::as<arma::cube>(eta_d[c]);
      arma::mat Q = create_Q(eta_c.n_rows + 1);
      gamma(d, c) = eta_to_gamma(eta_c, Q);
    }
  }
  return gamma;
}
