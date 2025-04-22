//Viterbi algorithm for NHMM and MHMM, single_sequence
#ifndef VITERBI_NHMM_H
#define VITERBI_NHMM_H

#include "config.h"

double univariate_viterbi(
    arma::uvec& q, const arma::vec& log_pi, const arma::cube& log_A, 
    const arma::mat& log_py
);

#endif
