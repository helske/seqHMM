#' Forward and Backward Probabilities for Hidden Markov Model
#'
#' The `forward_backward` function computes forward and backward 
#' probabilities of a hidden Markov model.
#'
#' @export
#' @param model A hidden Markov model.
#' @param forward_only If `TRUE`, only forward probabilities are computed. The 
#' default is `FALSE`.
#' @param log_space If `TRUE` (default), forward and backward probabilities are 
#' computed and returned in logarithmic scale for numerical stability. If 
#' `FALSE`, scaling is used instead, which is somewhat faster. Not used for 
#' `nhmm` and `mnhmm` objects which always use logarithmic scale.
#' @param threads Number of threads used in parallel computing. The default 
#' is `1`. Not used for `nhmm` and `mnhmm` objects.
#' @param as_data_frame If `TRUE` (default), the output is returned as a 
#' data.frame. Otherwise, a list of array(s) is returned. Ignored if 
#' `log_space` is `FALSE`, in which case list of arrays is always returned.
#' @param ... Ignored.
#' @return If `as_data_frame` is `TRUE` a `data.frame` with 
#' log-values of forward and backward probabilities. If `FALSE` or 
#' `log_space = FALSE`, a list with components
#' * forward_probs\cr If `log_space = FALSE`, scaled forward probabilities, 
#'   i.e. probability of state given observations up to that time point. 
#'   If `log_space = TRUE`, logarithms of non-scaled forward probabilities.
#' * backward_probs\cr Scaled backward probabilities (`log_space = FALSE`),
#'   or logarithms of non-scaled backward probabilities(`log_space = TRUE`).
#' * scaling_factors\cr Sum of non-scaled forward probabilities at each time 
#'   point. Only computed if `log_space = FALSE`.
#' In case of multiple observations, these are computed independently for each 
#' sequence.
#' @examples
#' # Load a pre-defined MHMM
#' data("mhmm_biofam")
#'
#' # Compute forward and backward probabilities
#' fb <- forward_backward(mhmm_biofam, as_data_frame = FALSE)
#'
#' # The most probable hidden state at time t
#' # given the observations up to time t for the first subject:
#' apply(fb$forward_probs[, , 1], 2, which.max)
#' 
#' @rdname forward_backward
#' @export
forward_backward <- function(model, ...) {
  UseMethod("forward_backward", model)
}
#' @rdname forward_backward
#' @export
forward_backward.hmm <- function(model, forward_only = FALSE, 
                                 log_space = TRUE, threads = 1,
                                 as_data_frame = TRUE, ...) {
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)
  if (!log_space) {
    out <- forwardbackward(
      model$transition_probs, emissionArray,
      model$initial_probs, obsArray, forward_only, threads
    )
  } else {
    out <- log_forwardbackward(
      model$transition_probs, emissionArray,
      model$initial_probs, obsArray, forward_only, threads
    )
  }
  if (is.null(time_names <- colnames(model$observations[[1]]))) {
    time_names <- seq_len(model$length_of_sequences)
  }
  if (is.null(sequence_names <- rownames(model$observations[[1]]))) {
    sequence_names <- seq_len(model$n_sequences)
  }
  dimnames(out$forward_probs) <- list(
    "state" = model$state_names,
    "time" = time_names, "sequence" = sequence_names
  )
  if (!log_space) {
    dimnames(out$scaling_factors) <- list("time" = time_names, "sequence" = sequence_names)
  }
  if (!forward_only) {
    dimnames(out$backward_probs) <- list(
      "state" = model$state_names,
      "time" = time_names, "sequence" = sequence_names
    )
  }
  if (as_data_frame && log_space) {
    do.call(rbind, lapply(names(out), function(x) {
      cbind(expand.grid(dimnames(out[[x]])), log_probability = c(out[[x]]),
            variable = x)
    }))
  } else {
    out
  }
}
#' @rdname forward_backward
#' @export
forward_backward.mhmm <- function(model, forward_only = FALSE, 
                                  log_space = TRUE, threads = 1,
                                  as_data_frame = TRUE, ...) {
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  model <- .combine_models(model)
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)
  if (!log_space) {
    out <- forwardbackwardx(
      model$transition_probs, emissionArray,
      model$initial_probs, obsArray, model$coefficients, model$X, model$n_states_in_clusters, forward_only, threads
    )
  } else {
    out <- log_forwardbackwardx(
      model$transition_probs, emissionArray,
      model$initial_probs, obsArray, model$coefficients, model$X, model$n_states_in_clusters, forward_only, threads
    )
  }
  if (is.null(time_names <- colnames(model$observations[[1]]))) {
    time_names <- seq_len(model$length_of_sequences)
  }
  if (is.null(sequence_names <- rownames(model$observations[[1]]))) {
    sequence_names <- seq_len(model$n_sequences)
  }
  dimnames(out$forward_probs) <- list(
    "state" = model$state_names,
    "time" = time_names, "sequence" = sequence_names
  )
  if (!log_space) {
    dimnames(out$scaling_factors) <- list("time" = time_names, "sequence" = sequence_names)
  }
  if (!forward_only) {
    dimnames(out$backward_probs) <- list(
      "state" = model$state_names,
      "time" = time_names, "sequence" = sequence_names
    )
  }
  if (as_data_frame && log_space) {
    do.call(rbind, lapply(names(out), function(x) {
      cbind(expand.grid(dimnames(out[[x]])), log_probability = c(out[[x]]),
            variable = x)
    }))
  } else {
    out
  }
}
#' @rdname forward_backward
#' @export
forward_backward.nhmm <- function(model, forward_only = FALSE, 
                                  as_data_frame = TRUE, ...) {
  X_initial <- t(model$X_initial)
  X_transition <- aperm(model$X_transition, c(3, 1, 2))
  X_emission <- aperm(model$X_emission, c(3, 1, 2))
  obsArray <- create_obsArray(model)
  beta_i_raw <- stan_to_cpp_initial(
    model$estimation_results$parameters$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    model$estimation_results$parameters$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    model$estimation_results$parameters$beta_o_raw,
    1,
    model$n_channels > 1
  )
  out <- list()
  if (model$n_channels == 1) {
    out$forward_probs <- forward_nhmm_singlechannel(
      beta_i_raw, X_initial,
      beta_s_raw, X_transition,
      beta_o_raw, X_emission,
      array(obsArray[1, , ], dim(obsArray)[2:3]))
    if (!forward_only) {
      out$backward_probs <- backward_nhmm_singlechannel(
        beta_i_raw, X_initial,
        beta_s_raw, X_transition,
        beta_o_raw, X_emission,
        array(obsArray[1, , ], dim(obsArray)[2:3]))
    }
    if (is.null(time_names <- colnames(model$observations))) {
      time_names <- seq_len(model$length_of_sequences)
    }
    if (is.null(sequence_names <- rownames(model$observations))) {
      sequence_names <- seq_len(model$n_sequences)
    }
  } else {
    out$forward_probs <- forward_nhmm_multichannel(
      beta_i_raw, X_initial,
      beta_s_raw, X_transition,
      beta_o_raw, X_emission,
      obsArray, model$n_symbols)
    if (!forward_only) {
      out$backward_probs <- backward_nhmm_multichannel(
        beta_i_raw, X_initial,
        beta_s_raw, X_transition,
        beta_o_raw, X_emission,
        obsArray, model$n_symbols)
    }
    if (is.null(time_names <- colnames(model$observations[[1]]))) {
      time_names <- seq_len(model$length_of_sequences)
    }
    if (is.null(sequence_names <- rownames(model$observations[[1]]))) {
      sequence_names <- seq_len(model$n_sequences)
    }
  }
  dimnames(out$forward_probs) <- list(
    "state" = model$state_names,
    "time" = time_names, "sequence" = sequence_names
  )
  if (!forward_only) {
    dimnames(out$backward_probs) <- list(
      "state" = model$state_names,
      "time" = time_names, "sequence" = sequence_names
    )
  }
  if (as_data_frame) {
    do.call(rbind, lapply(names(out), function(x) {
      cbind(expand.grid(dimnames(out[[x]])), log_probability = c(out[[x]]),
            variable = x)
    }))
  } else {
    out
  }
}

#' @rdname forward_backward
#' @export
forward_backward.mnhmm <- function(model, forward_only = FALSE, 
                                   as_data_frame = TRUE, ...) {
  X_initial <- t(model$X_initial)
  X_transition <- aperm(model$X_transition, c(3, 1, 2))
  X_emission <- aperm(model$X_emission, c(3, 1, 2))
  X_cluster <- t(model$X_cluster)
  obsArray <- create_obsArray(model)
  S <- model$n_states
  D <- model$n_clusters
  beta_i_raw <- stan_to_cpp_initial(
    model$estimation_results$parameters$beta_i_raw,
    model$n_clusters
  )
  beta_s_raw <- stan_to_cpp_transition(
    model$estimation_results$parameters$beta_s_raw,
    model$n_clusters
  )
  beta_o_raw <- stan_to_cpp_emission(
    model$estimation_results$parameters$beta_o_raw,
    model$n_clusters,
    model$n_channels > 1
  )
  theta_raw <- model$estimation_results$parameters$theta_raw
  out <- list()
  if (model$n_channels == 1) {
    out$forward_probs <- forward_mnhmm_singlechannel(
      beta_i_raw, X_initial,
      beta_s_raw, X_transition,
      beta_o_raw, X_emission,
      theta_raw, X_cluster,
      array(obsArray, dim(obsArray)[2:3]))
    if (!forward_only) {
      out$backward_probs <- backward_mnhmm_singlechannel(
        beta_i_raw, X_initial,
        beta_s_raw, X_transition,
        beta_o_raw, X_emission,
        theta_raw, X_cluster,
        array(obsArray, dim(obsArray)[2:3]))
    }
    if (is.null(time_names <- colnames(model$observations))) {
      time_names <- seq_len(model$length_of_sequences)
    }
    if (is.null(sequence_names <- rownames(model$observations))) {
      sequence_names <- seq_len(model$n_sequences)
    }
  } else {
    out$forward_probs <- forward_mnhmm_multichannel(
      beta_i_raw, X_initial,
      beta_s_raw, X_transition,
      beta_o_raw, X_emission,
      theta_raw, X_cluster,
      obsArray, model$n_symbols)
    if (!forward_only) {
      out$backward_probs <- backward_mnhmm_multichannel(
        beta_i_raw, X_initial,
        beta_s_raw, X_transition,
        beta_o_raw, X_emission,
        theta_raw, X_cluster,
        obsArray, model$n_symbols)
    }
    if (is.null(time_names <- colnames(model$observations[[1]]))) {
      time_names <- seq_len(model$length_of_sequences)
    }
    if (is.null(sequence_names <- rownames(model$observations[[1]]))) {
      sequence_names <- seq_len(model$n_sequences)
    }
  }
  state_names <- paste0(
    rep(model$cluster_names, each = model$n_states), ": ",
    model$state_names
  )
  dimnames(out$forward_probs) <- list(
    "state" = state_names, "time" = time_names, "id" = sequence_names
  )
  if (!forward_only) {
    dimnames(out$backward_probs) <- list(
      "state" = state_names, "time" = time_names, "id" = sequence_names
    )
  }
  if (as_data_frame) {
    do.call(rbind, lapply(names(out), function(x) {
      cbind(expand.grid(dimnames(out[[x]])), log_probability = c(out[[x]]),
            variable = x)
    }))
  } else {
    out
  }
}
