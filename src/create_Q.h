#ifndef CREATE_Q_H
#define CREATE_Q_H

#include "config.h"

arma::mat create_Q(const arma::uword n);
arma::field<arma::mat> create_Q(const arma::uvec n);

#endif
