#' Simulate Trajectories with FAN-HMM
#'
#' Simulate new sequences of observed and hidden states given an estimated FAN-HMM
#' and potentially new covariate data.
#'
#' @param model An object of class `fanhmm`.
#' @param newdata A data frame used to update `model`.
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
  W <- update_W_for_fanhmm(model, model$data)
  if (identical(coefs, "random")) {
    coefs <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL)
  } else {
    if (is.null(coefs$initial_probs)) coefs$initial_probs <- NULL
    if (is.null(coefs$transition_probs)) coefs$transition_probs <- NULL
    if (is.null(coefs$emission_probs)) coefs$emission_probs <- NULL
  }
  model$etas <- stats::setNames(
    create_initial_values(coefs, model, init_sd), c("pi", "A", "B")
  )
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  model$gammas$B <- eta_to_gamma_cube(model$etas$B)
  out <- simulate_fanhmm_singlechannel(
    model$etas$pi, model$X_pi, model$etas$A, W$W_A, model$etas$B, W$W_B,
    as.integer(obs_1) - 1, !is.null(autoregression_formula)
  )
  for (i in seq_len(model$n_sequences)) {
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
  out$observations[] <- symbol_names[c(out$observations) + 1]
  model$observations <- suppressWarnings(suppressMessages(
    seqdef(t(out$observations), alphabet = symbol_names, cnames = seq_len(T_))
  ))
  model$data[[response_name]] <- factor(
    c(t(model$observations)), 
    levels = TraMineR::alphabet(model$observations)
  )
  if (!is.null(autoregression_formula)) {
    model$data[[paste0("lag_", response_name)]] <- group_lag(model$data, id, response_name)
  }
  attr(model$X_pi, "X_mean") <- TRUE
  attr(model$X_A, "X_mean") <- TRUE
  attr(model$X_B, "X_mean") <- TRUE
  attr(model$X_pi, "R_inv") <- NULL
  attr(model$X_A, "R_inv") <- NULL
  attr(model$X_B, "R_inv") <- NULL
  model <- update(model, model$data)
  tQs <- t(create_Q(n_states))
  
  if (!attr(model$X_pi, "icpt_only")) {
    coef_names <- attr(model$X_pi, "coef_names")
    model$gammas$pi <- gamma_to_gamma_std(
      model$gammas$pi, solve(attr(model$X_pi, "R_inv")), 
      coef_names, attr(model$X_pi, "X_mean")
    )
    model$etas$pi[] <- tQs %*% model$gammas$pi
  }
  if (!attr(model$X_A, "icpt_only")) {
    coef_names <- attr(model$X_A, "coef_names")
    model$gammas$A <- gamma_to_gamma_std(
      model$gammas$A, solve(attr(model$X_A, "R_inv")), 
      coef_names, attr(model$X_A, "X_mean")
    )
    for (s in seq_len(n_states)) {
      model$etas$A[, , s] <- tQs %*% model$gammas$A[, , s]
    }
  }
  if (!attr(model$X_B, "icpt_only")) {
    tQm <- t(create_Q(model$n_symbols))
    coef_names <- attr(model$X_B, "coef_names")
    model$gammas$B <- gamma_to_gamma_std(
      model$gammas$B, solve(attr(model$X_B, "R_inv")), 
      coef_names, attr(model$X_B, "X_mean")
    )
    for (s in seq_len(n_states)) {
      model$etas$B[, , s] <- tQm %*% model$gammas$B[, , s]
    }
  }
  list(model = model, states = states)
}
