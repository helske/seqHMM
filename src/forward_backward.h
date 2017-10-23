#ifndef FB_H
#define FB_H

#include <RcppArmadillo.h>

void uvForward(const arma::sp_mat& transition_t, const arma::cube& emission, const arma::vec& init,
  const arma::umat& obs, arma::mat& alpha, arma::vec& scales);
void uvBackward(const arma::sp_mat& transition, const arma::cube& emission,
  const arma::umat& obs, arma::mat& beta, const arma::vec& scales);

void internalForwardx(const arma::sp_mat& transition_t, const arma::cube& emission,
                      const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, arma::mat& scales, unsigned int threads);
void internalBackwardx(const arma::sp_mat& transition, const arma::cube& emission,
                       const arma::ucube& obs, arma::cube& beta, const arma::mat& scales, unsigned int threads);

void internalForward(const arma::mat& transition, const arma::cube& emission,
                     const arma::vec& init, const arma::ucube& obs, arma::cube& alpha, arma::mat& scales, unsigned int threads);
void internalBackward(const arma::mat& transition, const arma::cube& emission,
                      const arma::ucube& obs, arma::cube& beta, const arma::mat& scales, unsigned int threads);

#endif
