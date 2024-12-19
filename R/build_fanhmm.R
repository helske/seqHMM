#' Build a Feedback-Augmented Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_fanhmm <- function(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, autoregression_formula, 
    feedback_formula, data, time, id, state_names = NULL) {
  
  stopifnot(
    !is.na(autoregression_formula) || !is.na(feedback_formula),
    "Provide {.arg autoregression_formula} and/or {.arg feedback_formula} for FAN-HMM."
  )
  out <- create_base_nhmm(
    observations, data, time, id, n_states, state_names, channel_names = NULL,
    initial_formula, transition_formula, emission_formula) 
  out[c("cluster_names", "n_clusters", "X_omega")] <- NULL
  setNames(
    create_initial_values(list(), out$model, 0), c("pi", "A", "B")
  )
  stopifnot_(
    out$n_channels == 1L,
    "Currently only single-channel responses are supported for FAN-HMM.")
  if(is.na(autoregression_formula)) {
    out$model$W_A <- array(
      0, c(0L, out$model$length_of_sequences, out$model$n_sequences)
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
    out$model$W_A <- W_A$X
    out$model$rhos$A <- create_rho_A_inits(
      NULL, n_states, out$model$n_symbols, nrow(out$model$W_A), 0
    )
    out$extras$intercept_only <- FALSE
    np_rho_A <- W_A$n_pars
  }
  if(is.na(feedback_formula)) {
    out$model$W_B <- array(
      0, c(0L, out$model$length_of_sequences, out$model$n_sequences)
    )
    out$model$rhos$B <- create_rho_B_inits(
      NULL, n_states, out$model$n_symbols, 0, 0
    )
    np_rho_B <- 0
  } else {
    W_B <- model_matrix_feedback_formula_formula(
      feedback_formula_formula, data, 
      out$model$n_sequences, 
      out$model$length_of_sequences, n_states,
      out$model$n_symbols, time, id, 
      out$model$sequence_lengths
    )
    out$model$W_B <- W_B$X
    out$model$rhos$B <- create_rho_B_inits(
      NULL, n_states, out$model$n_symbols, nrow(out$model$W_B), 0
    )
    out$extras$intercept_only <- FALSE
    np_rho_B <- W_B$n_pars
  }
  structure(
    c(
      out$model,
      list(autoregression_formula = autoregression_formula, 
           feedback_formula_formula = feedback_formula_formula
      )
    ),
    class = "fanhmm",
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
