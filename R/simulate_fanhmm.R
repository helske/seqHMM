#' Simulate Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' non-homogeneous hidden Markov model.
#'
#' @param n_sequences The number of sequences to simulate.
#' @param sequence_lengths The lengths of the simulated sequences. 
#' Either a scalar or vector of length `n_sequences`.  
#' @param n_symbols A scalar or vector of length `n_channels` giving the number 
#' of observed symbols per channel.
#' @param n_states The number of hidden states.
#' @param coefs If `coefs = "random"` (default), random coefficient values are 
#' used. Otherwise `coefs` should be named list of `eta_pi`, `eta_A`, 
#' and `eta_B`.
#' @param init_sd Standard deviation of the normal distribution used to generate
#' random coefficient values. Default is `2`.
#' @inheritParams estimate_mnhmm
#'
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_fanhmm <- function(
    n_sequences, sequence_lengths, n_symbols, n_states, 
    initial_formula = ~1, transition_formula = ~1, 
    emission_formula = ~1, autoregression_formula = ~1, 
    feedback_formula = ~1, obs_0, data, time, id, coefs = "random", 
    response_name = "y", init_sd = 2) {
  
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
    "{.arg data}  is missing, use {.fn simulate_mhmm} instead."
  )
  sequence_lengths <- rep(sequence_lengths, length.out = n_sequences)
  n_channels <- length(n_symbols)
  stopifnot_(
    n_channels == 1,
    "Currently only single-channel responses are supported for FAN-HMM."
  )
  stopifnot_(
    !missing(obs_0) && 
      checkmate::test_integerish(x = obs_0, lower = 1L, upper = n_symbols),
    "Argument {.arg obs_0} should be an integer vector of length {.val n_sequences}."
  )
  symbol_names <- as.character(seq_len(n_symbols))
  T_ <- max(sequence_lengths)
  data[[response_name]] <- factor(rep(symbol_names, length = nrow(data)), 
                                  levels = symbol_names)
  data0 <- data[data[[time]] == min(data[[time]]), ]
  data0[[time]] <- min(data[[time]]) - min(diff(sort(unique(data[[time]]))))
  data0[[response_name]] <- obs_0
  data0 <- rbind(
    data0,
    data
  )
 
  model <- build_fanhmm(
    response_name, n_states, initial_formula, transition_formula, emission_formula, 
    autoregression_formula, feedback_formula, data0, time, id, scale = FALSE
  )
  
  X_A <- X_B <- vector("list", n_symbols)
  for (i in seq_len(n_symbols)) {
    data[[response_name]] <- factor(symbol_names[i], levels = symbol_names)
    data[[paste0("lag_", response_name)]] <- data[[response_name]]
    mod <- update(model, data)
    X_A[[i]] <- mod$X_A
    X_B[[i]] <- mod$X_B
  }
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
  model$etas <- setNames(
    create_initial_values(coefs, model, init_sd), c("pi", "A", "B")
  )
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  model$gammas$B <- eta_to_gamma_cube(model$etas$B)
  out <- simulate_fanhmm_singlechannel(
    model$etas$pi, model$X_pi, 
    model$etas$A, X_A, 
    model$etas$B, X_B,
    as.integer(obs_0) - 1
  )
  for (i in seq_len(model$n_sequences)) {
    Ti <- sequence_lengths[i]
    if (Ti < T_) {
      out$states[(Ti + 1):T_, i] <- NA
      out$observations[, (Ti + 1):T_, i] <- NA
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
  list(model = model, states = states)
}
