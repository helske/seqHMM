#' Simulate Trajectories with FAN-HMM
#'
#' Simulate new sequences of observed and hidden states given an estimated FAN-HMM
#' and potentially new covariate data.
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
  obs_1 <- model$observations[, 1]
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
  state_names <- model$state_names
  out$states[] <- state_names[c(out$states) + 1]
  states <- suppressWarnings(suppressMessages(
    seqdef(
      matrix(
        t(out$states),
        n_sequences, max(sequence_lengths)
      ), 
      alphabet = state_names, cnames = seq_len(T_)
    )
  ))
  out$observations[] <- model$symbol_names[c(out$observations) + 1]
  model$observations <- suppressWarnings(suppressMessages(
    seqdef(t(out$observations), alphabet = model$symbol_names, cnames = seq_len(T_))
  ))
  response_name <- model$channel_names[1]
  model$data[[response_name]] <- factor(
    c(t(model$observations)), 
    levels = TraMineR::alphabet(model$observations)
  )
  model <- update(model, model$data)
  list(model = model, states = states)
}
