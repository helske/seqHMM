#ifndef LIST_TO_FIELD_H
#define LIST_TO_FIELD_H

#include "config.h"

// Converts a D x C nested lists of 3D arrays into a 2D field of cubes
// for eta_B and gamma_B in MNHMM case
arma::field<arma::cube> cubelist_to_2d_field(const Rcpp::List& x);

// Converts a nested lists of matrices into a 2D field of matrices
// for X_B (C lists of N lists, channels and IDs)
// for W_A (P lists of N lists, combinations and IDs)
arma::field<arma::mat> matlist_to_2d_field(const Rcpp::List& x);

// Converts a P x C x N nested list of matrices into a 3D field of matrices
// for W_B (P combinations, C channels, N individuals)
arma::field<arma::mat> matlist_to_3d_field(const Rcpp::List& x);

// Converts a P x C x N nested list of vectors into a 3D field of vectors
// for W_X_B (P combinations, C channels, N individuals)
arma::field<arma::vec> veclist_to_3d_field(const Rcpp::List& x);

#endif
