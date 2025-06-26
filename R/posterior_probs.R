#' Posterior Probabilities for Hidden Markov Models
#'
#' Function `posterior_probs` computes the posterior probabilities of hidden 
#' states of a (mixture) hidden Markov model.
#' @export
#' @param model A hidden Markov model object.
#' @param ... Ignored.
#' @return A data frame of posterior probabilities for each state and sequence.
#' @examples
#' # Load a pre-defined MHMM
#' data("mhmm_biofam")
#'
#' # Compute posterior probabilities
#' pb <- posterior_probs(mhmm_biofam)
#'
#' @rdname posterior_probs
#' @export
posterior_probs <- function(model, ...) {
  UseMethod("posterior_probs", model)
}
#' @rdname posterior_probs
#' @export
posterior_probs.hmm <- function(model, ...) {
  # avoid CRAN check warning due to NSE
  time <- id <- ll <- probability <- log_alpha <- log_beta <- NULL
  out <- forward_backward(model)
  out[, ll := logSumExp(log_alpha[time == time[.N]]), by = id, 
      showProgress = FALSE]
  out[, probability := exp(log_alpha + log_beta - ll)][, -c("log_alpha", "log_beta", "ll")]
}
#' @rdname posterior_probs
#' @export
posterior_probs.mhmm <- function(model, ...) {
  # avoid CRAN check warning due to NSE
  time <- id <- ll <- probability <- log_alpha <- log_beta <- NULL
  out <- forward_backward(model)
  out[, ll := logSumExp(log_alpha[time == time[.N]]), by = id, 
      showProgress = FALSE]
  out[, probability := exp(log_alpha + log_beta - ll)][, -c("log_alpha", "log_beta", "ll")]
}
#' @rdname posterior_probs
#' @export
posterior_probs.nhmm <- function(model, ...) { 
  # avoid CRAN check warning due to NSE
  ll <- probability <- log_alpha <- log_beta <- NULL
  out <- forward_backward(model)
  time_var <- model$time_variable
  id_var <- model$id_variable
  out[, ll := logSumExp(log_alpha[time_var == time_var[.N]]), 
      by = id_var, env = list(id_var = id_var, time_var = time_var), 
      showProgress = FALSE]
  out[, probability := exp(log_alpha + log_beta - ll)]
  out[, -c("log_alpha", "log_beta", "ll")]
}
#' @rdname posterior_probs
#' @export
posterior_probs.mnhmm <- function(model, ...) {
  # avoid CRAN check warning due to NSE
  ll <- probability <- log_alpha <- log_beta <- NULL
  out <- forward_backward(model)
  time_var <- model$time_variable
  id_var <- model$id_variable
  out[, ll := logSumExp(log_alpha[time_var == time_var[.N]]), 
      by = id_var, env = list(id_var = id_var, time_var = time_var), 
      showProgress = FALSE]
  out[, probability := exp(log_alpha + log_beta - ll)]
  out[, -c("log_alpha", "log_beta", "ll")]
}
