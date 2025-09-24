#ifndef SOFTMAX_H
#define SOFTMAX_H

#include "config.h"

arma::vec softmax(const arma::vec& x);
arma::vec log_softmax(const arma::vec& x);

#endif
