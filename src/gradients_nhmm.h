#ifndef GRADIENTS_NHMM_H
#define GRADIENTS_NHMM_H

#include "config.h"

void gradient_wrt_pi(
    arma::mat& grad, arma::mat& tmpmat,
    const arma::mat& log_py, const arma::mat& log_beta, const double ll, 
    const arma::vec& pi, const arma::mat& X, const arma::uword i);
void gradient_wrt_A(
    arma::mat& grad, arma::mat& tmpmat,
    const arma::mat& log_py, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& A, 
    const arma::cube& X, const arma::uword i, const arma::uword t, 
    const arma::uword s);
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::ucube& obs, const arma::vec& log_pi, const arma::mat& log_beta, 
    const double ll, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::field<arma::cube>& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c);
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::ucube& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll,  const arma::cube& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::field<arma::cube>& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c);

#endif
