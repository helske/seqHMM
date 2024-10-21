#' Number of Observations in Hidden Markov Model
#'
#' Extract the number of non-missing observations of HMM. When computing nobs 
#' for a multichannel model with $C$ channels, each observed value in a single 
#' channel amounts to $1/C$ observation, i.e. a fully observed time point for 
#' a single sequence amounts to one observation.
#' @param object An object of class `hmm`, `mhmm`, `nhmm`, or `mnhmm`.
#' @param ... Ignored.
#' @rdname nobs
#' @export
nobs.hmm <- function(object, ...) {
  attr(object, "nobs")
}
#' @rdname nobs
#' @export
nobs.mhmm <- function(object, ...) {
  attr(object, "nobs")
}
#' @rdname nobs
#' @export
nobs.nhmm <- function(object, ...) {
  attr(object, "nobs")
}
#' @rdname nobs
#' @export
nobs.mnhmm <- function(object, ...) {
  attr(object, "nobs")
}