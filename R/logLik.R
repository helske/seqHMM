#' Log-likelihood of the Hidden Markov Model
#'
#'
#' @param object Ahidden Markov model.
#' @param partials Return a vector containing the individual contributions of 
#' each sequence to the total log-likelihood. The default is `FALSE`, which 
#' returns the sum of all log-likelihood components.
#' @param threads Number of threads to use in parallel computing. The default 
#' is 1.
#' @param log_space Make computations using log-space instead of scaling for 
#' greater numerical stability at the cost of decreased computational 
#' performance. The default is `TRUE`.
#' @param ... Ignored.
#' @return Log-likelihood of the hidden Markov model. This is an object of class
#' `logLik` with attributes `nobs` and `df` inherited from the model object.
#' @rdname logLik
#' @export
logLik.hmm <- function(object, partials = FALSE, threads = 1, 
                       log_space = TRUE, ...) {
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  df <- attr(object, "df")
  nobs <- attr(object, "nobs")
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
  structure(
    if (partials) ll else sum(ll), 
    class = "logLik", 
    df = df, nobs = nobs
  )
}
#' @rdname logLik
#' @export
logLik.mhmm <- function(object, partials = FALSE, threads = 1, 
                        log_space = TRUE, ...) {
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  df <- attr(object, "df")
  nobs <- attr(object, "nobs")
  object <- .combine_models(object)
  obsArray <- create_obsArray(object)
  emissionArray <- create_emissionArray(object)
  if (!log_space) {
    ll <- logLikMixHMM(
      object$transition_probs, emissionArray, object$initial_probs, obsArray,
      object$coefficients, object$X, object$n_states_in_clusters, threads
    )
  } else {
    ll <- log_logLikMixHMM(
      object$transition_probs, emissionArray, object$initial_probs, obsArray,
      object$coefficients, object$X, object$n_states_in_clusters, threads
    )
  }
  structure(
    if (partials) ll else sum(ll), 
    class = "logLik", 
    df = df, nobs = nobs
  )
}
#' @rdname logLik
#' @export
logLik.nhmm <- function(object, partials = FALSE, ...) {
  df <- attr(object, "df")
  nobs <- attr(object, "nobs")
  out <- forward_backward(object, forward_only = TRUE)
  ll <- apply(out$forward_probs[, object$length_of_sequences, ], 2, logSumExp)
  structure(
    if (partials) ll else sum(ll), 
    class = "logLik", 
    df = df, nobs = nobs
  )
}
#' @rdname logLik
#' @export
logLik.mnhmm <- function(object, partials = FALSE, ...) {
  df <- attr(object, "df")
  nobs <- attr(object, "nobs")
  out <- forward_backward(object, forward_only = TRUE)
  ll <- apply(out$forward_probs[, object$length_of_sequences, ], 2, logSumExp)
  structure(
    if (partials) ll else sum(ll), 
    class = "logLik", 
    df = df, nobs = nobs
  )
}