#' The seqHMM package
#'
#' The seqHMM package is designed for fitting hidden (or latent) Markov models (HMMs) and
#' mixture hidden Markov models (MHMMs) for social sequence data and other categorical
#' time series. The package supports models for one or multiple subjects with one or
#' multiple interdependent sequences (channels). External covariates can be added to
#' explain cluster membership in mixture models. The package provides functions for evaluating
#' and comparing models, as well as functions for easy plotting of multichannel sequences
#' and hidden Markov models. Common restricted versions of (M)HMMs are also supported,
#' namely Markov models, mixture Markov models, and latent class models.
#'
#' Maximum likelihood estimation via the EM algorithm and direct numerical maximization
#' with analytical gradients is supported. All main algorithms are written in C++.
#' Parallel computation is implemented via OpenMP for pre-2.0.0 functions, while 
#' estimation of non-homogenous models support parallelization via `future` 
#' package by parallelization of multistart optimizations and bootstrap sampling.
#'
#' @docType package
#' @name seqHMM-package
#' @aliases seqHMM
#' @useDynLib seqHMM, .registration = TRUE
#' @import igraph
#' @import gridBase
#' @import grid
#' @import nloptr
#' @import Rcpp
#' @importFrom Rcpp loadModule evalCpp
#' @importFrom stats logLik vcov predict update coef nobs
#' @importFrom utils head modifyList
#' @importFrom TraMineR seqlegend
#' @importFrom graphics par plot.new
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_ribbon geom_line facet_wrap
#' @importFrom collapse fndistinct seq_row rsplit qM
#' @import data.table
#' @references Helske S. and Helske J. (2019). Mixture Hidden Markov Models for Sequence Data: The seqHMM Package in R,
#' Journal of Statistical Software, 88(3), 1-32. doi:10.18637/jss.v088.i03
"_PACKAGE"
#' @importFrom TraMineR alphabet
#' @export
TraMineR::alphabet
#' @importFrom TraMineR seqstatf
#' @export
TraMineR::seqstatf
#' @importFrom TraMineR seqdef
#' @export
TraMineR::seqdef
