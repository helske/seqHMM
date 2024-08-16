#ifndef GETPROBS_H
#define GETPROBS_H

#include <RcppArmadillo.h>

arma::mat get_pi(
    const arma::mat& beta_raw, const arma::mat& X, const int logspace
);
arma::field<arma::cube> get_A(
    const arma::cube& beta_raw, const arma::cube& X, const int logspace
);
arma::field<arma::cube> get_B(
    const arma::cube& beta_raw, const arma::cube& X, const int logspace
);
arma::field<arma::cube> get_multichannel_B(
    const arma::vec& beta_raw, const arma::cube& X, unsigned int S, 
    unsigned int C, const arma::uvec& M, const int logspace
);

#endif