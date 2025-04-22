#include "list_to_field.h"

// Converts a D x C nested lists of 3D arrays into a 2D field of cubes
// for eta_B and gamma_B in MNHMM case
arma::field<arma::cube> cubelist_to_2d_field(const Rcpp::List& x) {
  arma::uword D = x.size();
  Rcpp::List list = x[0];
  arma::uword C = list.size();
  arma::field<arma::cube> out(D, C);
  for (arma::uword d = 0; d < D; ++d) {
    list = x[d];
    for (arma::uword c = 0; c < C; ++c) {
      out(d, c) = Rcpp::as<arma::cube>(list[c]);
    }
  }
  return out;
}
// Converts a nested lists of matrices into a 2D field of matrices
// for X_B (C lists of N lists, channels and IDs)
// for W_A (P lists of N lists, combinations and IDs)
arma::field<arma::mat> matlist_to_2d_field(const Rcpp::List& x) {
  arma::uword C = x.size();
  Rcpp::List list = x[0];
  arma::uword N = list.size();
  arma::field<arma::mat> out(C, N);
  for (arma::uword c = 0; c < C; ++c) {
    list = x[c];
    for (arma::uword n = 0; n < N; ++n) {
      out(c, n) = Rcpp::as<arma::mat>(list[n]);
    }
  }
  return out;
}

// Converts a P x C x N nested list of matrices into a 3D field of matrices
// for W_B (P combinations, C channels, N individuals)
arma::field<arma::mat> matlist_to_3d_field(const Rcpp::List& x) {
  arma::uword P = x.size();
  Rcpp::List list1 = x[0];
  arma::uword C = list1.size();
  Rcpp::List list2 = list1[0];
  arma::uword N = list2.size();
  arma::field<arma::mat> out(P, C, N);
  for (arma::uword p = 0; p < P; ++p) {
    for (arma::uword c = 0; c < C; ++c) {
      list1 = x[p];
      list2 = list1[c];
      for (arma::uword n = 0; n < N; ++n) {
        out(p, c, n) = Rcpp::as<arma::mat>(list2[n]);
      }
    }
  }
  return out;
}

// Converts a P x C x N nested list of vectors into a 3D field of vectors
// for W_X_B (P combinations, C channels, N individuals)
arma::field<arma::vec> veclist_to_3d_field(const Rcpp::List& x) {
  arma::uword P = x.size();
  Rcpp::List list1 = x[0];
  arma::uword C = list1.size();
  Rcpp::List list2 = list1[0];
  arma::uword N = list2.size();
  arma::field<arma::vec> out(P, C, N);
  for (arma::uword p = 0; p < P; ++p) {
    for (arma::uword c = 0; c < C; ++c) {
      list1 = x[p];
      list2 = list1[c];
      for (arma::uword n = 0; n < N; ++n) {
        out(p, c, n) = Rcpp::as<arma::vec>(list2[n]);
      }
    }
  }
  return out;
}