#ifndef SEQHMM_H
#define SEQHMM_H

// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::vec reparma(const arma::vec& x, const arma::uvec& y);

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

unsigned int optCoef(arma::mat& weights, const arma::ucube& obs, const arma::cube& emission,
                     const arma::mat& bsi, arma::mat& coef, const arma::mat& X,
                     const arma::uvec& cumsumstate, const arma::uvec& numberOfStates, int trace);

arma::vec gCoef(const arma::ucube& obs, const arma::mat& bsi, const arma::cube& emission,
                const arma::mat& weights, const arma::mat& X, const arma::uvec& cumsumstate,
                const arma::uvec& numberOfStates);

arma::mat hCoef(const arma::mat& weights, const arma::mat& X);

//log-space versions
double logSumExp(const arma::vec& x);

void log_internalForwardx(const arma::mat& transition, const arma::cube& emission,
                          const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, unsigned int threads);

void log_internalForward(const arma::mat& transition, const arma::cube& emission,
                         const arma::vec& init, const arma::ucube& obs, arma::cube& alpha, unsigned int threads);

void log_internalBackward(const arma::mat& transition, const arma::cube& emission,
                          const arma::ucube& obs, arma::cube& beta, unsigned int threads);

unsigned int log_optCoef(arma::mat& weights, const arma::ucube& obs, const arma::cube& emission, const arma::mat& initk,
                         const arma::cube& beta, const arma::vec& ll, arma::mat& coef, const arma::mat& X,
                         const arma::uvec& cumsumstate, const arma::uvec& numberOfStates, int trace);

arma::vec log_gCoef(const arma::ucube& obs, const arma::cube& beta, const arma::cube& emission, const arma::mat& initk,
                    const arma::mat& weights, const arma::vec& ll, const arma::mat& X, const arma::uvec& cumsumstate,
                    const arma::uvec& numberOfStates);

void log_internalBackward_single(const arma::mat& transition, const arma::cube& emission,
                                 const arma::ucube& obs, arma::mat& beta, int k);

void log_internalForwardx_single(const arma::mat& transition, const arma::cube& emission,
                                 const arma::vec& init, const arma::ucube& obs, arma::mat& alpha, int k);

#endif
