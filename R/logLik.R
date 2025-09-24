#' Log-likelihood of a Hidden Markov Model
#'
#' @param object A hidden Markov model.
#' @param partials Return a vector containing the individual contributions of 
#' each sequence to the total log-likelihood. The default is `FALSE`, which 
#' returns the sum of all log-likelihood components.
#' @param threads Number of threads to use in parallel computing. The default 
#' is 1.
#' @param log_space Make computations using log-space instead of scaling for 
#' greater numerical stability at the cost of decreased computational 
#' performance. The default is `FALSE`.
#' @param ... Ignored.
#' @return Log-likelihood of the hidden Markov model. This is an object of class
#' `logLik` with attributes `nobs` and `df` inherited from the model object.
#' @rdname logLik_hmm
#' @export
logLik.hmm <- function(object, partials = FALSE, threads = 1, 
                       log_space = FALSE, ...) {
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
#' @rdname logLik_hmm
#' @export
logLik.mhmm <- function(object, partials = FALSE, threads = 1, 
                        log_space = FALSE, ...) {
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
#' Log-likelihood of a Non-homogeneous Hidden Markov Model
#'
#' @param object A hidden Markov model.
#' @param partials Return a vector containing the individual contributions of 
#' each sequence to the total log-likelihood. The default is `FALSE`, which 
#' returns the sum of all log-likelihood components.
#' @param ... Ignored.
#' @return Log-likelihood of the hidden Markov model. This is an 
#' object of class `logLik` with attributes `nobs` and `df` inherited from the 
#' model object.
#' @rdname logLik_nhmm
#' @export
#' @export
logLik.nhmm <- function(object, partials = FALSE, ...) {
  # Avoid warnings due to NSE
  ll <- time <- id <- log_alpha <- NULL
  df <- attr(object, "df")
  nobs <- attr(object, "nobs")
  if (partials || is.null(object$estimation_results$loglik)) {
    if (inherits(object, "fanhmm")) {
      ll <- Rcpp_loglik_fanhmm(
        create_obs(object), object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, 
        io(object$X_pi), io(object$X_A), io(object$X_B),
        iv(object$X_A), iv(object$X_B),
        tv(object$X_A), tv(object$X_B),
        object$etas$eta_pi, object$etas$eta_A, object$etas$eta_B,
        object$prior_obs, object$W_X_B
      )
    } else {
      ll <- Rcpp_loglik_nhmm(
        create_obs(object), object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, 
        io(object$X_pi), io(object$X_A), io(object$X_B),
        iv(object$X_A), iv(object$X_B),
        tv(object$X_A), tv(object$X_B),
        object$etas$eta_pi, object$etas$eta_A, object$etas$eta_B
      )
    }
    if (!is.null(object$estimation_results$lambda)) {
      ll <- ll - 0.5 * object$estimation_results$lambda * sum(unlist(object$etas)^2) / object$n_sequences
    }
  } else {
    ll <- object$estimation_results$loglik
  }
  structure(
    if (partials) ll else sum(ll), 
    class = "logLik", 
    df = df, nobs = nobs
  )
}
#' @rdname logLik_nhmm
#' @export
logLik.mnhmm <- function(object, partials = FALSE, ...) {
  # Avoid warnings due to NSE
  ll <- time <- id <- log_alpha <- NULL
  df <- attr(object, "df")
  nobs <- attr(object, "nobs")
  if (partials || is.null(object$estimation_results$loglik)) {
    if (inherits(object, "fanhmm")) {
      ll <- Rcpp_loglik_mfanhmm(
        create_obs(object), object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, object$X_omega,
        io(object$X_pi), io(object$X_A), io(object$X_B), io(object$X_omega),
        iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
        object$etas$eta_pi, object$etas$eta_A, object$etas$eta_B, 
        object$etas$eta_omega, object$prior_obs, object$W_X_B
      )
    } else {
      ll <- Rcpp_loglik_mnhmm(
        create_obs(object), object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, object$X_omega,
        io(object$X_pi), io(object$X_A), io(object$X_B), io(object$X_omega),
        iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
        object$etas$eta_pi, object$etas$eta_A, object$etas$eta_B, 
        object$etas$eta_omega
      )
    }
    if (!is.null(object$estimation_results$lambda)) {
      ll <- ll - 0.5 * object$estimation_results$lambda * sum(unlist(object$etas)^2) / object$n_sequences
    }
  } else {
    ll <- object$estimation_results$loglik
  }
  structure(
    if (partials) ll else sum(ll), 
    class = "logLik", 
    df = df, nobs = nobs
  )
}
