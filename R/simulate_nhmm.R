#' Simulate Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' non-homogeneous hidden Markov model.
#' 
#' @param n_states The number of hidden states.
#' @param data A data frame containing the variables used in the model 
#' formulas. Note that this should also include also the response variable(s), 
#' which are used to define the number of observed symbols (using [levels()]) 
#' and the length of sequences. The actual values of the response variables 
#' does not matter though, as they are replaced by the simulated values. The 
#' exception is the first time point in FAN-HMM case: If the `emission_formula` 
#' contains lagged responses, the response variable values at the first time 
#' point are used to define the emissions at the second time point, and the 
#' simulations are done from the second time point onwards. This matches the 
#' case `prior_obs = "fixed"` in [estimate_nhmm()].
#' @param coefs Same as argument `inits` in [estimate_nhmm()]. If `NULL` 
#' (default), the model parameters are generated randomly. If you want to 
#' simulate new sequences based on an estimated model `fit`, you can use 
#' `coefs = fit$etas` and `init_sd = 0`.
#' @param init_sd Standard deviation of the normal distribution used to generate
#' random coefficient values. Default is `2`. Setting this to zero allows you 
#' to use the coefficients fixed in `coefs`.
#' @inheritParams estimate_nhmm
#'
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_nhmm <- function(
    n_states, emission_formula, initial_formula = ~1, transition_formula = ~1, 
    data, id, time, coefs = NULL, init_sd = 2) {
  
  stopifnot_(
    !missing(data),
    "{.arg data} is missing, use {.fn simulate_hmm} instead."
  )
  stopifnot_(
    !missing(emission_formula),
    "Argument {.arg emission_formula} is missing."
  )
  if (inherits(emission_formula, "formula")) {
    responses <- get_responses(emission_formula)
    C <- length(responses)
    if (C > 1L) {
      rhs <- deparse1(emission_formula[[3L]])
      emission_formula <- lapply(
        responses, \(y) stats::as.formula(
          paste(y, " ~ ", rhs), 
          env = environment(emission_formula)
        )
      )
    } else {
      emission_formula <- list(emission_formula)
    }
  } else {
    responses <- vapply(emission_formula, get_responses, allow_mv = FALSE, "")
    C <- length(responses)
  }
  data <- .check_data(data, id, time, responses)
  for (y in responses) {
    l <- levels(data[[y]])
    data[, y := ifelse(is.na(y[1]), l[1], y[1]), by = id, env = list(y = y)]
  }
  if (is.null(coefs)) {
    coefs <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL)
  } else {
    if (is.null(coefs$initial_probs)) coefs$initial_probs <- NULL
    if (is.null(coefs$transition_probs)) coefs$transition_probs <- NULL
    if (is.null(coefs$emission_probs)) coefs$emission_probs <- NULL
  }
  model <- build_nhmm(
    n_states, emission_formula, initial_formula, transition_formula, 
    data, id, time, coefs = coefs
  )
  if (inherits(model, "fanhmm")) {
    W <- update_W_for_fanhmm(model)
    out <- Rcpp_simulate_fanhmm(
      create_obs(model), model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$prior_obs, model$W_X_B, W$W_A, W$W_B
    )
  } else {
    out <- Rcpp_simulate_nhmm(
      create_obs(model), model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B
    )
  }
  for (i in seq_len(model$n_channels)) {
    y <- unlist(lapply(out$observations, \(y) y[i, ] + 1L))
    model$data[[model$responses[i]]][] <- 
      model$symbol_names[[i]][y]
  }
  states <- cbind(
    model$data[, list(id, time), env = list(id = id, time = time)], 
    state = model$state_names[unlist(out$states) + 1L]
  )
  list(model = model, states = states)
}
