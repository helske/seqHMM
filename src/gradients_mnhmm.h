#ifndef GRADIENTS_MNHMM_H
#define GRADIENTS_MNHMM_H

#include "config.h"

void gradient_wrt_omega(
    arma::mat& grad, arma::mat& tmpmat,
     const arma::vec& omega, 
    const arma::vec& loglik_i, const arma::vec& loglik, const arma::mat& X, 
    const arma::uword i);
void gradient_wrt_pi(
    arma::mat& grad, arma::mat& tmpmat,
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::vec>& pi, const arma::mat& X, const arma::uword i,
    const arma::uword d);
void gradient_wrt_A(
    arma::mat& grad, arma::mat& tmpmat,
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_alpha, const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube> A,
    const arma::cube& X, const arma::uword i, const arma::uword t, 
    const arma::uword s, const arma::uword d);
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec,
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::field<arma::vec>& log_pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::field<arma::cube>& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c, const arma::uword d);
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::field<arma::cube>& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c,
    const arma::uword d);

#endif
