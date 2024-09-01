#ifndef BACKWARD_NHMM_H
#define BACKWARD_NHMM_H

#include <RcppArmadillo.h>

arma::mat univariate_backward_nhmm(
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py);

#endif