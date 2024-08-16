#' Log-likelihood of the Hidden Markov Model
#'
#' Function `logLik.hmm` computes the log-likelihood value of a hidden Markov model.
#'
#'
#' @export
#' @param object A  hidden Markov model of class `hmm`.
#' @param partials Return a vector containing the individual contributions of each sequence to the total log-likelihood.
#'   The default is `FALSE`, which returns the sum of all log-likelihood components.
#' @param threads Number of threads to use in parallel computing. The default is 1.
#' @param log_space Make computations using log-space instead of scaling for greater
#' numerical stability at the cost of decreased computational performance.
#'   The default is `TRUE`.
#' @param ... Ignored.
#' @return Log-likelihood of the hidden Markov model. This is an object of class
#' `logLik` with attributes `nobs` and `df` inherited from the model object.
#' @seealso [build_hmm()] and [fit_model()] for building and
#'   fitting Hidden Markov models.
logLik.hmm <- function(object, partials = FALSE, threads = 1, log_space = FALSE, ...) {
  check_positive_integer(threads)
  obsArray <- create_obsArray(object)
  emissionArray <- create_emissionArray(object)

  if (!log_space) {
    ll <- logLikHMM(
      object$transition_probs, emissionArray,
      object$initial_probs, obsArray, threads
    )
  } else {
    ll <- log_logLikHMM(
      object$transition_probs, emissionArray,
      object$initial_probs, obsArray, threads
    )
  }

  structure(if (partials) ll else sum(ll), class = "logLik", df = attr(object, "df"), nobs = attr(object, "nobs"))
}
