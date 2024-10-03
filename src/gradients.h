#ifndef GRADIENTS_H
#define GRADIENTS_H

#include <RcppArmadillo.h>
arma::mat gradient_wrt_pi(
    const arma::mat& log_py, const arma::mat& log_beta, const double ll, 
    const arma::vec& Pi, const arma::mat& X, const unsigned int i);

arma::mat gradient_wrt_pi(
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::vec>& Pi, const arma::mat& X, const unsigned int i,
    const unsigned int d);

arma::mat gradient_wrt_A(
    const arma::mat& log_py, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& A, 
    const arma::cube& X, const unsigned int i, const unsigned int t, 
    const unsigned int s);

arma::mat gradient_wrt_A(
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_alpha, const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube> A,
    const arma::cube& X, const unsigned int i, const unsigned int t, 
    const unsigned int s, const unsigned int d);

// NHMM SC
arma::mat gradient_wrt_B_t0(
    const arma::umat& obs, const arma::vec& log_Pi, const arma::mat& log_beta, 
    const double ll, const arma::cube& B, const arma::cube& X, 
    const unsigned int i, const unsigned int s);
arma::mat gradient_wrt_B(
    const arma::umat& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, 
    const arma::cube& log_A, const arma::cube& B, const arma::cube& X, 
    const unsigned int i, const unsigned int s, const unsigned int t);
// NHMM MC
arma::mat gradient_wrt_B_t0(
    const arma::ucube& obs, const arma::vec& log_Pi, const arma::mat& log_beta, 
    const double ll, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const unsigned int i, const unsigned int s, 
    const unsigned int c);
arma::mat gradient_wrt_B(
    const arma::ucube& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll,  const arma::cube& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const unsigned int i, 
    const unsigned int s, const unsigned int t, const unsigned int c);
// MNHMM SC
arma::mat gradient_wrt_B_t0(
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::field<arma::vec>& log_Pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& B, const arma::cube& X, 
    const unsigned int i, const unsigned int s, const unsigned d);
arma::mat gradient_wrt_B(
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, const arma::field<arma::cube>& B, 
    const arma::cube& X, const unsigned int i, const unsigned int s, 
    const unsigned int t, const unsigned int d);
// MNHMM MC
arma::mat gradient_wrt_B_t0(
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::field<arma::vec>& log_Pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const unsigned int i, const unsigned int s, 
    const unsigned int c, const unsigned int d);
arma::mat gradient_wrt_B(
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const unsigned int i, 
    const unsigned int s, const unsigned int t, const unsigned int c,
    const unsigned int d);

#endif
