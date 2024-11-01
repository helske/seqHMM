#ifndef CREATEQ_H
#define CREATEQ_H

#include <RcppArmadillo.h>
arma::vec2 givens(const double a, const double b);
arma::mat compute_cs(const arma::uword n);
arma::mat create_Q(const arma::uword n);
#endif
