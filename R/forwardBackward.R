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
  d[, c("cluster", "state") := tstrsplit(state, ":", fixed = TRUE)]
  d
}
#' @rdname forward_backward
#' @export
forward_backward.nhmm <- function(model, forward_only = FALSE,  ...) {
  obs <- create_obsArray(model)
  fp <- forward_nhmm(
    obs, model$sequence_lengths, model$n_symbols, 
    model$X_pi, model$X_A, model$X_B, 
    io(model$X_pi), io(model$X_A), io(model$X_B),
    iv(model$X_A), iv(model$X_B),
    tv(model$X_A), tv(model$X_B),
    model$etas$pi, model$etas$A, model$etas$B
  )
  if (!forward_only) {
    bp <- backward_nhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B),
      tv(model$X_A), tv(model$X_B),
      model$etas$pi, model$etas$A, model$etas$B
    )
  }
  ids <- unique(model$data[[model$id_variable]])
  times <- unique(model$data[[model$time_variable]])
  d <- data.table(
    expand.grid(
      state = model$state_names,
      time = times,
      id = ids,
      stringsAsFactors = FALSE
    )[, 3:1],
    log_alpha = c(fp),
    log_beta = if (forward_only) NULL else c(bp),
    key = c("id", "time")
  )
  setnames(d, c("id", "time"), c(model$id_variable, model$time_variable))
  stats::na.omit(d)
}

#' @rdname forward_backward
#' @export
forward_backward.mnhmm <- function(model, forward_only = FALSE,  ...) {
  # avoid CRAN check warnings due to NSE
  state <- NULL
  obs <- create_obsArray(model)
  fp <- forward_mnhmm(
    obs, model$sequence_lengths, model$n_symbols, 
    model$X_pi, model$X_A, model$X_B, model$X_omega,
    io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
    iv(model$X_A), iv(model$X_B),
    tv(model$X_A), tv(model$X_B),
    model$etas$pi, model$etas$A, model$etas$B, model$etas$omega
  )
  if (!forward_only) {
    bp <- backward_mnhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B),
      tv(model$X_A), tv(model$X_B),
      model$etas$pi, model$etas$A, model$etas$B, model$etas$omega
    )
  }
  state_names <- paste0(
    rep(model$cluster_names, each = model$n_states), ": ",
    unlist(model$state_names)
  )
  ids <- unique(model$data[[model$id_variable]])
  times <- unique(model$data[[model$time_variable]])
  
  d <- data.table(
    expand.grid(
      state = state_names,
      time = times,
      id = ids,
      stringsAsFactors = FALSE
    )[, 3:1],
    log_alpha = c(fp),
    log_beta = if (forward_only) NULL else c(bp),
    key = c("id", "time")
  )
  setnames(d, c("id", "time"), c(model$id_variable, model$time_variable))
  d[, c("cluster", "state") := tstrsplit(state, ":", fixed = TRUE)]
  stats::na.omit(d)
}
