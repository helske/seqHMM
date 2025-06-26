#' Forward and Backward Probabilities for Hidden Markov Model
#'
#' The `forward_backward` function computes forward and backward 
#' probabilities of a hidden Markov model.
#'
#' @export
#' @param model A hidden Markov model.
#' @param forward_only If `TRUE`, only forward probabilities are computed. The 
#' default is `FALSE`.
#' @param ... Ignored.
#' @return A `data.frame` with log-values of forward and backward probabilities. 
#' @examples
#' # Load a pre-defined MHMM
#' data("mhmm_biofam")
#'
#' # Compute forward and backward probabilities
#' fb <- forward_backward(mhmm_biofam)
#'
#' head(fb)
#' 
#' @rdname forward_backward
#' @export
forward_backward <- function(model, ...) {
  UseMethod("forward_backward", model)
}
#' @rdname forward_backward
#' @export
forward_backward.hmm <- function(model, forward_only = FALSE, ...) {
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)
  out <- log_forwardbackward(
    model$transition_probs, emissionArray,
    model$initial_probs, obsArray, forward_only, 1L
  )
  if (is.null(times <- colnames(model$observations[[1]]))) {
    times <- seq_len(model$length_of_sequences)
  }
  if (is.null(ids <- rownames(model$observations[[1]]))) {
    ids <- seq_len(model$n_sequences)
  }
  ids <- as_factor(ids)
  d <- data.table(
    expand.grid(
      state = model$state_names,
      time = times,
      id = as.factor(ids),
      stringsAsFactors = FALSE
    )[, 3:1],
    log_alpha = c(out$forward_probs),
    log_beta = if (forward_only) NULL else c(out$backward_probs),
    key = c("id", "time")
  )
  d
}
#' @rdname forward_backward
#' @export
forward_backward.mhmm <- function(model, forward_only = FALSE, ...) {
  # avoid CRAN check warnings due to NSE
  state <- NULL
  model <- .combine_models(model)
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)
  out <- log_forwardbackwardx(
    model$transition_probs, emissionArray,
    model$initial_probs, obsArray, model$coefficients, model$X, 
    model$n_states_in_clusters, forward_only, 1L
  )
  if (is.null(times <- colnames(model$observations[[1]]))) {
    times <- seq_len(model$length_of_sequences)
  }
  if (is.null(ids <- rownames(model$observations[[1]]))) {
    ids <- seq_len(model$n_sequences)
  }
  ids <- as_factor(ids)
  d <- data.table(
    expand.grid(
      state = model$state_names,
      time = times,
      id = ids,
      stringsAsFactors = FALSE
    )[, 3:1],
    log_alpha = c(out$forward_probs),
    log_beta = if (forward_only) NULL else c(out$backward_probs),
    key = c("id", "time")
  )
  d[, c("cluster", "state") := tstrsplit(state, ": ", fixed = TRUE)]
  setcolorder(d, c("id", "time", "cluster", 
                   setdiff(names(d), c("id", "time", "cluster"))))
}
#' @rdname forward_backward
#' @export
forward_backward.nhmm <- function(model, forward_only = FALSE,  ...) {
  obs <- create_obs(model)
  if (inherits(model, "fanhmm")) {
    fp <- Rcpp_forward_fanhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B),
      tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B,
      model$prior_obs, model$W_X_B
    )
    if (!forward_only) {
      bp <- Rcpp_backward_fanhmm(
        obs, model$sequence_lengths, model$n_symbols, 
        model$X_pi, model$X_A, model$X_B, 
        io(model$X_pi), io(model$X_A), io(model$X_B),
        iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
        model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B,
        model$prior_obs, model$W_X_B
      )
    }
  } else {
    fp <- Rcpp_forward_nhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B),
      tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B
    )
    if (!forward_only) {
      bp <- Rcpp_backward_nhmm(
        obs, model$sequence_lengths, model$n_symbols, 
        model$X_pi, model$X_A, model$X_B, 
        io(model$X_pi), io(model$X_A), io(model$X_B),
        iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
        model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B
      )
    }
  }
  S <- model$n_states
  id <- model$id_variable
  time <- model$time_variable
  d <- model$data[rep(seq_row(model$data), each = S), 
                  list(id, time), 
                  env = list(id = id, time = time, S = S)]
  n <- nrow(d)
  set(d, j = "state", value = rep_len(model$state_names, n))
  set(d, j = "log_alpha", value = unlist(fp))
  if (!forward_only) {
    set(d, j = "log_beta", value = unlist(bp))
  }
  d[]
}

#' @rdname forward_backward
#' @export
forward_backward.mnhmm <- function(model, forward_only = FALSE,  ...) {
  obs <- create_obs(model)
  if (inherits(model, "fanhmm")) {
    fp <- Rcpp_forward_mfanhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$gammas$gamma_omega,
      model$prior_obs, model$W_X_B
    )
    if (!forward_only) {
      bp <- Rcpp_backward_mfanhmm(
        obs, model$sequence_lengths, model$n_symbols, 
        model$X_pi, model$X_A, model$X_B, model$X_omega,
        io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
        iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
        model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
        model$gammas$gamma_omega,
        model$prior_obs, model$W_X_B
      )
    }
  } else {
    fp <- Rcpp_forward_mnhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$gammas$gamma_omega
    )
    if (!forward_only) {
      bp <- Rcpp_backward_mnhmm(
        obs, model$sequence_lengths, model$n_symbols, 
        model$X_pi, model$X_A, model$X_B, model$X_omega,
        io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
        iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
        model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
        model$gammas$gamma_omega
      )
    }
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
  set(d, j = "log_alpha", value = unlist(fp))
  if (!forward_only) {
    set(d, j = "log_beta", value = unlist(bp))
  }
  d[]
}
