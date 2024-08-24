#ifndef FORWARD_NHMM_H
#define FORWARD_NHMM_H

#include <RcppArmadillo.h>

double univariate_forward_nhmm(
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py, 
    arma::subview_col<unsigned int> q);

#endif