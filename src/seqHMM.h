#ifndef SEQHMM_H
#define SEQHMM_H


#ifdef _OPENMP
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#endif

#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#endif
