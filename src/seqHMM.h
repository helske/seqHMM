#ifndef SEQHMM_H
#define SEQHMM_H

//#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

void internalForward2(const arma::mat& transition, const arma::cube& emission, 
  const arma::vec& init, const arma::Cube<int>& obs, arma::cube& alpha);

double logSumExp(const arma::vec& x);
arma::vec reparma(arma::vec x, Rcpp::IntegerVector y);

void internalForward(const arma::mat& transition, const arma::cube& emission, 
  const arma::vec& init, const arma::Cube<int>& obs, arma::cube& alpha);
void internalForwardx(const arma::mat& transition, const arma::cube& emission,
  const arma::mat& init, const arma::Cube<int>& obs, arma::cube& alpha);
void internalBackward(const arma::mat& transition, const arma::cube& emission, 
  const arma::Cube<int>& obs, arma::cube& beta);


arma::mat optCoef(const arma::icube& obs, const arma::cube& emission, const arma::mat& initk, 
  const arma::cube& beta, const arma::vec& ll, arma::mat& coef, const arma::mat& X, 
  const Rcpp::IntegerVector cumsumstate, const Rcpp::IntegerVector numberOfStates, int trace);

arma::vec gCoef(const arma::icube& obs, const arma::cube& beta, const arma::cube& emission, const arma::mat& initk,
  const arma::mat& weights, const arma::vec& ll, const arma::mat& X, const Rcpp::IntegerVector cumsumstate, 
  const Rcpp::IntegerVector numberOfStates);

arma::mat hCoef(const arma::mat& weights, const arma::mat& X);

#endif
