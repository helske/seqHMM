#ifndef SEQHMM_H
#define SEQHMM_H

#ifdef _OPENMP
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#endif

#define ARMA_NO_DEBUG
#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#endif
