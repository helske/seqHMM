#ifndef SUMTOZERO_H
#define SUMTOZERO_H

#include "config.h"

arma::mat sum_to_zero(const arma::mat& x, const arma::mat& Q);
arma::vec sum_to_zero(const arma::vec& x, const arma::mat& Q);

#endif
