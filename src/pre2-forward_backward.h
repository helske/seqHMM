#ifndef FB_H
#define FB_H

#include "config.h"

void uvForward(
    const arma::mat& transition_t, const arma::cube& emission, 
    const arma::vec& init, const arma::umat& obs, arma::mat& alpha, 
    arma::vec& scales
);

void uvBackward(
    const arma::mat& transition, const arma::cube& emission,
    const arma::umat& obs, arma::mat& beta, const arma::vec& scales
);

void internalForward(
    const arma::mat& transition_t, const arma::cube& emission,
    const arma::vec& init, const arma::ucube& obs, arma::cube& alpha, 
    arma::mat& scales, arma::uword threads
);

void internalForward(
    const arma::mat& transition_t, const arma::cube& emission,
    const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, 
    arma::mat& scales, arma::uword threads
);

void internalBackward(
    const arma::mat& transition, const arma::cube& emission,
    const arma::ucube& obs, arma::cube& beta, const arma::mat& scales, 
    arma::uword threads
);

#endif
