#ifndef VITERBI_NHMM_H
#define VITERBI_NHMM_H

#include <RcppArmadillo.h>

double univariate_viterbi_nhmm(
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py, 
    arma::subview_col<unsigned int> q);

#endif