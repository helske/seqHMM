#' Simulate Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' non-homogeneous hidden Markov model.
#'
#' @param n_sequences The number of sequences to simulate.
#' @param sequence_lengths The lengths of the simulated sequences. 
#' Either a scalar or vector of length `n_sequences`.  
#' @param n_symbols A scalar or a vector giving the number of observed symbols 
#' per response variable.
#' @param n_states The number of hidden states.
#' @param coefs If `NULL` (default), model parameters `eta_pi`, `eta_A`, and 
#' `eta_B` are generated randomly. You can also specify the coefficients by 
#' providing a named list with elements `eta_pi`, `eta_A`, and `eta_B`. 
#' Alternatively, values `initial_probs`, `transition_probs`, and
#' `emission_probs` can be used, which define the values of the \eqn{\eta} 
#' coefficients corresponding to the intercept terms, while rest of the values 
#' are generated randomly.
#' @param init_sd Standard deviation of the normal distribution used to generate
#' random coefficient values. Default is `2`. Setting this to zero allows you 
#' to use the coefficients fixed in `coefs`.
#' @inheritParams estimate_nhmm
#'
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_nhmm <- function(
    n_sequences, sequence_lengths, n_symbols, n_states, 
    emission_formula, initial_formula = ~1, transition_formula = ~1, 
    data, id, time, coefs = NULL, init_sd = 2) {
  
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
    "Argument {.arg n_symbols} must contain positive integer(s) larger than 1."
  )
  stopifnot_(
    !missing(n_states) && checkmate::test_int(x = n_states, lower = 2L), 
    "Argument {.arg n_states} must be a single positive integer larger than 1."
  )
  stopifnot_(
    !missing(data),
    "{.arg data}  is missing, use {.fn simulate_hmm} instead."
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
        responses, \(y) as.formula(
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
  stopifnot_(
    C == length(unique(responses)), 
    "Response names in {.arg responses} should be unique."
  )
  names(emission_formula) <- responses
  stopifnot_(
    all(idx <- !(responses %in% names(data))),
    "Variable{?s} with name{?s} {responses[!idx]} already found in {.arg data}.
    Use other name{?s} for the response variable{?s}."
  )
  stopifnot_(
    length(n_symbols) == length(responses),
    "`length(symbols)` is not equal to the number of responses 
    {length(responses)}."
  )
  sequence_lengths <- rep(sequence_lengths, length.out = n_sequences)
  n_channels <- length(n_symbols)
  T_ <- max(sequence_lengths)
  
  for (i in seq_len(n_channels)) {
    data[[responses[i]]] <- factor(1, levels = seq_len(n_symbols[i]))
  }
  model <- build_nhmm(
    n_states, emission_formula, initial_formula, transition_formula, 
    data, id, time, scale = TRUE
  )
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
  model$etas <- stats::setNames(
    create_initial_values(coefs, model, init_sd), c("pi", "A", "B")
  )
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
  out <- simulate_nhmm_cpp(
    model$etas$pi, model$X_pi, 
    model$etas$A, model$X_A, 
    model$etas$B, model$X_B, 
    model$n_symbols
  )
  for (i in seq_len(model$n_sequences)) {
    Ti <- sequence_lengths[i]
    if (Ti < T_) {
      out$states[(Ti + 1):T_, i] <- NA
      out$observations[, (Ti + 1):T_, i] <- NA
    }
  }
  for (i in seq_len(n_channels)) {
    model$data[[responses[i]]] <- factor(c(out$observations[i, , ]) + 1L)
  }
  states <- cbind(
    model$data[, list(id, time), env = list(id = id, time = time)], 
    state = model$state_names[c(out$states) + 1]
  )
  list(model = model, states = states)
}
