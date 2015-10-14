#ifndef SEQHMM_H
#define SEQHMM_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

void internalForward(const arma::mat& transition, const arma::cube& emission, 
  const arma::vec& init, const arma::icube& obs, arma::cube& alpha, arma::mat& scales);

arma::vec reparma(arma::vec x, Rcpp::IntegerVector y);

void internalForwardx(const arma::mat& transition, const arma::cube& emission,
  const arma::mat& init, const arma::icube& obs, arma::cube& alpha, arma::mat& scales);
void internalBackward(const arma::mat& transition, const arma::cube& emission, 
  const arma::icube& obs, arma::cube& beta, const arma::mat& scales);


unsigned int optCoef(arma::mat weights, const arma::icube& obs, const arma::cube& emission, const arma::mat& initk, 
  const arma::cube& beta, const arma::mat& scales, arma::mat& coef, const arma::mat& X, 
  const Rcpp::IntegerVector cumsumstate, const Rcpp::IntegerVector numberOfStates, int trace);

arma::vec gCoef(const arma::icube& obs, const arma::cube& beta, const arma::mat& scales, const arma::cube& emission, const arma::mat& initk,
  const arma::mat& weights, const arma::mat& X, const Rcpp::IntegerVector cumsumstate, 
  const Rcpp::IntegerVector numberOfStates);

arma::mat hCoef(const arma::mat& weights, const arma::mat& X);

#endif
