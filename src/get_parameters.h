#ifndef GETPROBS_H
#define GETPROBS_H

#include <RcppArmadillo.h>
#include "softmax.h"

arma::vec get_omega(const arma::mat& gamma_raw, const arma::vec& X);
arma::vec get_pi(const arma::mat& gamma_raw, const arma::vec& X);
arma::cube get_A(
    const arma::cube& gamma_raw, const arma::mat& X, const bool tv = true
  );
arma::cube get_B(
    const arma::cube& gamma_raw, const arma::mat& X,
    const bool add_missing = false, const bool tv = true
);
arma::field<arma::cube> get_B(
    const arma::field<arma::cube>& gamma_raw, const arma::mat& X, 
    const arma::uvec& M, const bool add_missing = false, 
    const bool tv = true
);

arma::vec get_log_omega(const arma::mat& gamma_raw, const arma::vec& X);
arma::vec get_log_pi(const arma::mat& gamma_raw, const arma::vec& X);
arma::cube get_log_A(
    const arma::cube& gamma_raw, const arma::mat& X, const bool tv = true
);
arma::cube get_log_B(
    const arma::cube& gamma_raw, const arma::mat& X,
    const bool add_missing = false, const bool tv = true
);
arma::field<arma::cube> get_log_B(
    const arma::field<arma::cube>& gamma_raw, const arma::mat& X, 
    const arma::uvec& M, const bool add_missing = false, 
    const bool tv = true
);
#endif
