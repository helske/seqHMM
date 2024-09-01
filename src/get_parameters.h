#ifndef GETPROBS_H
#define GETPROBS_H

#include <RcppArmadillo.h>

arma::mat get_omega(
    const arma::mat& theta_raw, const arma::mat& X, const int logspace
);
arma::mat get_pi(
    const arma::mat& beta_raw, const arma::mat& X, const int logspace
);
arma::field<arma::cube> get_A(
    const arma::cube& beta_raw, const arma::cube& X, const int logspace
);
arma::field<arma::cube> get_B(
    const arma::cube& beta_raw, const arma::cube& X, const int logspace, 
    const int add_missing = 0
);
arma::field<arma::cube> get_multichannel_B(
    const arma::vec& beta_raw, const arma::cube& X, unsigned int S, 
    unsigned int C, const arma::uvec& M, const int logspace, 
    const int add_missing = 0
);
arma::vec get_omega_i(
    const arma::mat& theta_raw, const arma::vec X, const int logspace
);
arma::vec get_pi_i(
    const arma::mat& beta_raw, const arma::vec X, const int logspace
);
arma::cube get_A_i(
    const arma::cube& beta_raw, const arma::mat& X, const int logspace
);
arma::cube get_B_i(
    const arma::cube& beta_raw, const arma::mat& X, const int logspace, 
    const int add_missing = 0
);
arma::field<arma::cube> get_multichannel_B_i(
    const arma::vec& beta_raw, const arma::mat& X, unsigned int S, 
    unsigned int C, const arma::uvec& M, const int logspace, 
    const int add_missing = 0
);
#endif