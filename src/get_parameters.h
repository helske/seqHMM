#ifndef GETPROBS_H
#define GETPROBS_H

#include <RcppArmadillo.h>

arma::vec get_omega(
    const arma::mat& gamma_raw, const arma::vec& X, const bool logspace
);
arma::vec get_pi(
    const arma::mat& gamma_raw, const arma::vec& X, const bool logspace
);
arma::cube get_A(
    const arma::cube& gamma_raw, const arma::mat& X, const bool logspace,
    const bool tv = true
);
arma::cube get_B(
    const arma::cube& gamma_raw, const arma::mat& X, const bool logspace, 
    const bool add_missing = false, const bool tv = true
);
arma::field<arma::cube> get_B(
    const arma::field<arma::cube>& gamma_raw, const arma::mat& X, 
    const arma::uvec& M, const bool logspace, const bool add_missing = false, 
    const bool tv = true
);
#endif
