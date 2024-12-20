#' Build a Feedback-Augmented Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_fanhmm <- function(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, autoregression_formula, 
    feedback_formula, data, time, id, state_names = NULL) {
  
  stopifnot_(
    !is.null(autoregression_formula) || !is.null(feedback_formula),
    "Provide {.arg autoregression_formula} and/or {.arg feedback_formula} for FAN-HMM."
  )
  stopifnot_(
    inherits(autoregression_formula, "formula"), 
    "Argument {.arg autoregression_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(feedback_formula, "formula"), 
    "Argument {.arg feedback_formula} must be a {.cls formula} object.")
  out <- create_base_nhmm(
    observations, data, time, id, n_states, state_names, channel_names = NULL,
    initial_formula, transition_formula, emission_formula)
  stopifnot_(
    !any(out$model$observations == attr(out$model$observations, "nr")),
    "FAN-HMM does not support missing values in the observations."
  )
  out[c("cluster_names", "n_clusters", "X_omega")] <- NULL
  out$model$etas <- setNames(
    create_initial_values(list(), out$model, 0), c("pi", "A", "B")
  )
  stopifnot_(
    out$model$n_channels == 1L,
    "Currently only single-channel responses are supported for FAN-HMM.")
  if(is.null(autoregression_formula)) {
    out$model$W_A <- array(
      0, c(0L, out$model$length_of_sequences - 1, out$model$n_sequences)
    )
    out$model$rhos$A <- create_rho_A_inits(
      NULL, n_states, out$model$n_symbols, 0, 0
    )
    np_rho_A <- 0
  } else {
    W_A <- model_matrix_autoregression_formula(
      autoregression_formula, data, 
      out$model$n_sequences, 
      out$model$length_of_sequences, n_states,
      out$model$n_symbols, time, id, 
      out$model$sequence_lengths
    )
    out$model$W_A <- W_A$X[, -1, , drop = FALSE]
    out$model$rhos$A <- create_rho_A_inits(
      NULL, n_states, out$model$n_symbols, nrow(out$model$W_A), 0
    )
    out$extras$intercept_only <- FALSE
    np_rho_A <- W_A$n_pars
  }
  if(is.null(feedback_formula)) {
    out$model$W_B <- array(
      0, c(0L, out$model$length_of_sequences - 1, out$model$n_sequences)
    )
    out$model$rhos$B <- create_rho_B_inits(
      NULL, n_states, out$model$n_symbols, 0, 0
    )
    np_rho_B <- 0
  } else {
    W_B <- model_matrix_feedback_formula(
      feedback_formula, data, 
      out$model$n_sequences, 
      out$model$length_of_sequences, n_states,
      out$model$n_symbols, time, id, 
      out$model$sequence_lengths
    )
    out$model$W_B <- W_B$X[, -1, , drop = FALSE]
    out$model$rhos$B <- create_rho_B_inits(
      NULL, n_states, out$model$n_symbols, nrow(out$model$W_B), 0
    )
    out$extras$intercept_only <- FALSE
    np_rho_B <- W_B$n_pars
  }
  out$model$obs_0 <- as.integer(out$model$observations[, 1]) - 1L
  out$model$observations <- out$model$observations[, -1]
  out$model$length_of_sequences <- out$model$length_of_sequences - 1
  out$model$sequence_lengths <- out$model$sequence_lengths - 1
  out$model$X_A <- out$model$X_A[, -1, , drop = FALSE]
  out$model$X_B <- out$model$X_B[, -1, , drop = FALSE]
  structure(
    c(
      out$model,
      list(
        autoregression_formula = autoregression_formula, 
        feedback_formula = feedback_formula
      )
    ),
    class = c("fanhmm", "nhmm"),
    nobs = attr(out$model$observations, "nobs"),
    df = out$extras$np_pi + out$extras$np_A + out$extras$np_B + np_rho_A + 
      np_rho_B,
    type = paste0(out$extras$multichannel, "fanhmm"),
    intercept_only = out$extras$intercept_only,
    np_pi = out$extras$np_pi,
    np_A = out$extras$np_A,
    np_B = out$extras$np_B,
    np_rho_A = np_rho_A,
    np_rho_B = np_rho_B
  )
}
