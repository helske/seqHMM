#ifndef SEQHMM_H
#define SEQHMM_H

#include <RcppArmadillo.h>

double logSumExp(const double& x,const double& y);

arma::vec reparma(arma::vec x, Rcpp::IntegerVector y);

void internalForwardMC(const arma::mat& transition, const arma::cube& emission, 
  const arma::vec& init, const arma::Cube<int>& obs, arma::cube& alpha);
void internalForward(const arma::mat& transition, const arma::mat& emission, 
  const arma::vec& init, const arma::Mat<int>& obs, arma::cube& alpha);
void internalForwardMCx(const arma::mat& transition, const arma::cube& emission,
  const arma::mat& init, const arma::Cube<int>& obs, arma::cube& alpha);
void internalForwardx(const arma::mat& transition, const arma::mat& emission,
  const arma::mat& init, const arma::imat& obs, arma::cube& alpha);

void internalBackwardMC(const arma::mat& transition, const arma::cube& emission, 
  const arma::Cube<int>& obs, arma::cube& beta);

void internalBackward(const arma::mat& transition, const arma::mat& emission, 
  const arma::Mat<int>& obs, arma::cube& beta);

void optCoef(arma::mat alpha, arma::mat beta, arma::vec& ll, arma::mat& coef,
  arma::mat& X, arma::mat& lweights,Rcpp::IntegerVector cumsumstate,
  Rcpp::IntegerVector numberOfStates);

arma::vec gCoef(arma::mat& usums, arma::mat& lweights, arma::mat& X);

arma::mat hCoef(arma::mat& usums, arma::mat& lweights, arma::mat& X);

void viterbiForEM(arma::mat& transition, arma::mat& emission, arma::vec& init, arma::imat& obs, 
  arma::vec& logp, arma::umat& q);

void viterbiForEMMC(arma::mat& transition, arma::cube& emission, arma::vec& init, arma::icube& obs, 
  arma::vec& logp, arma::umat& q);

void viterbiForEMx(arma::mat& transition, arma::mat& emission, arma::mat& init, arma::imat& obs, 
  arma::vec& logp, arma::umat& q);

void viterbiForEMMCx(arma::mat& transition, arma::cube& emission, arma::mat& init, arma::icube& obs, 
  arma::vec& logp, arma::umat& q);

#endif
