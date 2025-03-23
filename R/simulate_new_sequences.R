#' Simulate Trajectories with FAN-HMM
#'
#' Simulate new sequences of observed and hidden states given an estimated FAN-HMM
#' and potentially new covariate data, conditionally on the observations at 
#' the first time point.
#'
#' @param model An object of class `fanhmm`.
#' @param newdata Optional data frame used to update `model`.
#' @return A list with the model used in simulation with updated response variable, 
#' as well as the simulated hidden state sequences.
#' @export
simulate_new_sequences <- function(model, newdata = NULL) {
  
  stopifnot_(
    inherits(model, "fanhmm"),
    "Argument {.arg model} must an object of class {.cls fanhmm}."
  )
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  W <- update_W_for_fanhmm(model, newdata)
  times <- model$data[[model$time_variable]]
  obs_1 <- model$data[[model$response]][times == min(times)]
  out <- simulate_fanhmm_singlechannel(
    model$etas$pi, model$X_pi, model$etas$A, W$W_A, model$etas$B, W$W_B,
    as.integer(obs_1) - 1L, !is.null(model$autoregression_formula)
  )
  T_ <- model$length_of_sequences
  sequence_lengths <- model$sequence_lengths
  n_sequences <- model$n_sequences
  for (i in seq_len(n_sequences)) {
    Ti <- sequence_lengths[i]
    if (Ti < T_) {
      out$states[(Ti + 1):T_, i] <- NA
      out$observations[(Ti + 1):T_, i] <- NA
    }
  }
  model$data[[model$responses]] <- model$symbol_names[c(out$observations) + 1L]
  model <- update(model, model$data)
  states <- cbind(
    model$data[, c(model$id_variable, model$time_variable)],
    state = model$state_names[c(out$states) + 1L]
  )
  list(model = model, states = states)
}
