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
#' Parallel computation is implemented via OpenMP.
#' 
#' @docType package
#' @name seqHMM
#' @aliases seqHMM
#' @useDynLib seqHMM
#' @import igraph
#' @import gridBase
#' @import grid
#' @import nloptr
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix .bdiag
#' @importFrom stats logLik cmdscale complete.cases model.matrix BIC rnorm runif vcov
#' @importFrom TraMineR alphabet seqstatf seqdef seqlegend seqdist seqdistmc seqplot seqlength
#' @importFrom grDevices col2rgb rainbow
#' @importFrom graphics barplot par plot plot.new polygon strwidth text
#' @importFrom methods hasArg
#' @importFrom utils menu
NULL



