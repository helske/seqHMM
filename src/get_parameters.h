#ifndef GET_PROBS_H
#define GET_PROBS_H

#include "config.h"
#include "softmax.h"

arma::vec get_omega(
    const arma::mat& gamma, const arma::vec& X
);
arma::vec get_pi(
    const arma::mat& gamma, const arma::vec& X
);
arma::cube get_A(
    const arma::cube& gamma, const arma::mat& X, const bool tv = true
);
arma::cube get_B(
    const arma::cube& gamma, const arma::mat& X, const bool tv = true,
    const bool add_missing = false
);

#endif
