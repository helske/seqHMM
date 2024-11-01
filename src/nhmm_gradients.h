#ifndef GRADIENTS_H
#define GRADIENTS_H

#include <RcppArmadillo.h>

void gradient_wrt_omega(
    arma::mat& grad, arma::mat& tmpmat,
    const arma::mat& Q, const arma::vec& omega, 
    const arma::vec& loglik_i, const arma::vec& loglik, const arma::mat& X, 
    const arma::uword i);
  
void gradient_wrt_pi(
    arma::mat& grad, arma::mat& tmpmat,const arma::mat& Q,
    const arma::mat& log_py, const arma::mat& log_beta, const double ll, 
    const arma::vec& Pi, const arma::mat& X, const arma::uword i);

void gradient_wrt_pi(
    arma::mat& grad, arma::mat& tmpmat,const arma::mat& Q,
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::vec>& Pi, const arma::mat& X, const arma::uword i,
    const arma::uword d);

void gradient_wrt_A(
    arma::mat& grad, arma::mat& tmpmat,const arma::mat& Q,
    const arma::mat& log_py, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& A, 
    const arma::cube& X, const arma::uword i, const arma::uword t, 
    const arma::uword s);

void gradient_wrt_A(
    arma::mat& grad, arma::mat& tmpmat,const arma::mat& Q,
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_alpha, const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube> A,
    const arma::cube& X, const arma::uword i, const arma::uword t, 
    const arma::uword s, const arma::uword d);

// NHMM SC
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::umat& obs, const arma::vec& log_Pi, const arma::mat& log_beta, 
    const double ll, const arma::cube& B, const arma::cube& X, 
    const arma::uword i, const arma::uword s);
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::umat& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, 
    const arma::cube& log_A, const arma::cube& B, const arma::cube& X, 
    const arma::uword i, const arma::uword s, const arma::uword t);
// NHMM MC
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::ucube& obs, const arma::vec& log_Pi, const arma::mat& log_beta, 
    const double ll, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c);
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::ucube& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll,  const arma::cube& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c);
// MNHMM SC
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::field<arma::vec>& log_Pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uword i, const arma::uword s, const unsigned d);
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uword i, const arma::uword s, 
    const arma::uword t, const arma::uword d);
// MNHMM MC
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::field<arma::vec>& log_Pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c, const arma::uword d);
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, const arma::mat& Q,
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c,
    const arma::uword d);

#endif
