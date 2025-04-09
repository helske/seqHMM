#ifndef LIST_TO_FIELD_H
#define LIST_TO_FIELD_H

#include "config.h"

// Converts a list of lists of length D x C of arrays into a 2D field of cubes.
arma::field<arma::cube> list_to_2d_field(const Rcpp::List& x);

#endif