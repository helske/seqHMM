#' Simulate Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' non-homogeneous hidden Markov model.
#'
#' @param n_symbols A scalar giving the number of observed symbols.
#' @param autoregression_formula Formula for autoregression 
#' \eqn{y_t \to y_{t+1}}.
#' Default intercept means is shorthand for \eqn{y_{t+1} ~ y_{t}}, while 
#' additional terms in formula are interacted with the lagged responses. 
#' If `NULL`, no autoregression is used.
#' @param feedback_formula Formula for feedback \eqn{y_t \to z_{t+1}}.
#' Default intercept means is shorthand for \eqn{z_{t+1} ~ y_{t}}, while 
#' additional terms in formula are interacted with the lagged responses. 
#' If `NULL`, no feedback is used.
#' @param obs_1 Vector defining the observations at first time points in case 
#' of `autoregression_formula` is not `NULL`.
#' @inheritParams simulate_nhmm
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_fanhmm <- function(
    n_sequences, sequence_lengths, n_symbols, n_states, 
    initial_formula = ~1, transition_formula = ~1, 
    emission_formula = ~1, autoregression_formula = ~1, 
    feedback_formula = ~1, obs_1, data, id, time, responses, 
    coefs = NULL, init_sd = 2) {
  
  stopifnot_(
    !missing(n_sequences) && checkmate::test_int(x = n_sequences, lower = 1L), 
    "Argument {.arg n_sequences} must be a single positive integer."
  )
  stopifnot_(
    !missing(sequence_lengths) && 
      checkmate::test_integerish(x = sequence_lengths, lower = 2L), 
    "Argument {.arg sequence_lengths} must contain positive integer(s) larger than 1."
  )
  stopifnot_(
    !missing(n_symbols) && checkmate::test_integerish(x = n_symbols, lower = 2L), 
    "Argument {.arg n_symbols} must contain positive integer larger than 1."
  )
  stopifnot_(
    !missing(n_states) && checkmate::test_int(x = n_states, lower = 2L), 
    "Argument {.arg n_states} must be a single positive integer larger than 1."
  )
  stopifnot_(
    !missing(data),
    "{.arg data} is missing."
  )
  stopifnot_(
    checkmate::test_string(responses),
    "Argument {.arg responses} must be a scalar of type character."
  )
  stopifnot_(
    all(idx <- !(responses %in% names(data))),
    "Variable{?s} with name{?s} {responses[!idx]} found in {.arg data}.
    Use other name{?s} for the response variable{?s}."
  )
  sequence_lengths <- rep(sequence_lengths, length.out = n_sequences)
  n_channels <- length(responses)
  stopifnot_(
    n_channels == 1,
    "Currently only single-channel responses are supported for FAN-HMM."
  )
  if (!is.null(autoregression_formula)) {
    stopifnot_(
      !missing(obs_1) &&
        checkmate::test_integerish(x = obs_1, lower = 1L, upper = n_symbols),
      "Argument {.arg obs_1} should be an integer vector of length {.val n_sequences}."
    )
  }
  symbol_names <- as.character(seq_len(n_symbols))
  T_ <- max(sequence_lengths)
  data[[responses]] <- factor(rep(symbol_names, length = nrow(data)), 
                                  levels = symbol_names)
  
  model <- build_fanhmm(
    responses, n_states, initial_formula, transition_formula, emission_formula, 
    autoregression_formula, feedback_formula, data, id, time, scale = FALSE
  )
  
  W <- update_W_for_fanhmm(model)
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
  model$data[[responses]] <- factor(c(out$observations) + 1L)
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
  states <- cbind(
    model$data[, c(id, time)],
    state = model$state_names[c(out$states) + 1L]
  )
  list(model = model, states = states)
}
