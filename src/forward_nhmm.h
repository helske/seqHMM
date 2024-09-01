#ifndef FORWARD_NHMM_H
#define FORWARD_NHMM_H

#include <RcppArmadillo.h>

arma::mat univariate_forward_nhmm(
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py);

#endif