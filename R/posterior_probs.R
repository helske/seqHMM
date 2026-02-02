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
  obs <- create_obs(model)
  if (inherits(model, "fanhmm")) {
    pp <- Rcpp_posterior_probs_fanhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B),
      tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B,
      model$prior_obs, model$W_X_B
    )
  } else {
    pp <- Rcpp_posterior_probs_nhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B),
      tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B
    )
  }
  S <- model$n_states
  id <- model$id_variable
  time <- model$time_variable
  d <- model$data[rep(seq_row(model$data), each = S), 
                  list(id, time), 
                  env = list(id = id, time = time, S = S)]
  n <- nrow(d)
  set(d, j = "state", value = rep_len(model$state_names, n))
  set(d, j = "probability", value = unlist(pp))
  d[]
}
#' @rdname posterior_probs
#' @export
posterior_probs.mnhmm <- function(model, ...) {
  obs <- create_obs(model)
  if (inherits(model, "fanhmm")) {
    pp <- Rcpp_posterior_probs_mfanhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$gammas$gamma_omega
    )
  } else {
    pp <- Rcpp_posterior_probs_mnhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$gammas$gamma_omega
    )
  }

  D <- model$n_clusters
  S <- model$n_states
  DS <- D * S
  id <- model$id_variable
  time <- model$time_variable
  d <- model$data[rep(seq_row(model$data), each = DS), 
                  list(id, time), 
                  env = list(id = id, time = time, DS = DS)]
  n <- nrow(d)
  set(d, j = "cluster", value = rep_len(rep(model$cluster_names, each = S), n))
  set(d, j = "state", value = rep_len(unlist(model$state_names), n))
  set(d, j = "probability", value = unlist(pp))
  d[]
}
