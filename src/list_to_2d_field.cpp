#include "list_to_2d_field.h"
// Converts a list of lists of length D x C of arrays into a 2D field of cubes.
arma::field<arma::cube> list_to_2d_field(const Rcpp::List& x) {
  arma::uword D = x.size();
  Rcpp::List first_list = x[0];
  arma::uword C = first_list.size();
  arma::field<arma::cube> out(D, C);
  for (arma::uword d = 0; d < D; ++d) {
    Rcpp::List list = x[d];
    for (arma::uword c = 0; c < C; ++c) {
      out(d, c) = Rcpp::as<arma::cube>(list[c]);
    }
  }
  return out;
}