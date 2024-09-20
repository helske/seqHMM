#ifndef GETPROBS_H
#define GETPROBS_H

#include <RcppArmadillo.h>

arma::vec get_omega(
    const arma::mat& gamma_omega_raw, const arma::vec X, const int logspace
);
arma::vec get_pi(
    const arma::mat& beta_raw, const arma::vec X, const int logspace
);
arma::cube get_A(
    const arma::cube& beta_raw, const arma::mat& X, const int logspace
);
arma::cube get_B(
    const arma::cube& beta_raw, const arma::mat& X, const int logspace, 
    const int add_missing = 0
);
arma::field<arma::cube> get_B(
    const arma::field<arma::cube>& beta_raw, const arma::mat& X, 
    const arma::uvec& M, const int logspace, const int add_missing = 0
);
#endif
