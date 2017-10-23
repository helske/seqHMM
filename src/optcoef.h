#ifndef OPTCOEF_H
#define OPTCOEF_H


#include <RcppArmadillo.h>
unsigned int optCoef(arma::mat& weights, const arma::ucube& obs, const arma::cube& emission,
                     const arma::mat& bsi, arma::mat& coef, const arma::mat& X,
                     const arma::uvec& cumsumstate, const arma::uvec& numberOfStates, int trace);

arma::vec gCoef(const arma::ucube& obs, const arma::mat& bsi, const arma::cube& emission,
                const arma::mat& weights, const arma::mat& X, const arma::uvec& cumsumstate,
                const arma::uvec& numberOfStates);

arma::mat hCoef(const arma::mat& weights, const arma::mat& X);

//log-space versions

unsigned int log_optCoef(arma::mat& weights, const arma::ucube& obs, const arma::cube& emission, const arma::mat& initk,
                         const arma::cube& beta, const arma::vec& ll, arma::mat& coef, const arma::mat& X,
                         const arma::uvec& cumsumstate, const arma::uvec& numberOfStates, int trace);

arma::vec log_gCoef(const arma::ucube& obs, const arma::cube& beta, const arma::cube& emission, const arma::mat& initk,
                    const arma::mat& weights, const arma::vec& ll, const arma::mat& X, const arma::uvec& cumsumstate,
                    const arma::uvec& numberOfStates);


#endif
