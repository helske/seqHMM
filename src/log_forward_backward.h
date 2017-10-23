#ifndef LFB_H
#define LFB_H

#include <RcppArmadillo.h>

void log_internalForwardx(const arma::mat& transition, const arma::cube& emission,
                          const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, unsigned int threads);

void log_internalForward(const arma::mat& transition, const arma::cube& emission,
                         const arma::vec& init, const arma::ucube& obs, arma::cube& alpha, unsigned int threads);

void log_internalBackward(const arma::mat& transition, const arma::cube& emission,
                          const arma::ucube& obs, arma::cube& beta, unsigned int threads);


#endif
