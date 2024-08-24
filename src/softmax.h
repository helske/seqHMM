#ifndef SOFTMAX_H
#define SOFTMAX_H

#include <RcppArmadillo.h>
arma::vec softmax(const arma::vec& x, const int logspace);
#endif
