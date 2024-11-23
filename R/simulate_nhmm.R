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
simulate_nhmm <- function(
    n_sequences, sequence_lengths, n_symbols, n_states, 
    initial_formula = ~1, transition_formula = ~1, 
    emission_formula = ~1, data, time, id, coefs = "random", init_sd = 2) {
  
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
    "{.arg data}  is missing, use {.fn simulate_mhmm} instead."
  )
  sequence_lengths <- rep(sequence_lengths, length.out = n_sequences)
  n_channels <- length(n_symbols)
  symbol_names <- lapply(
    seq_len(n_channels), function(i) {
      as.character(seq_len(n_symbols[i]))
    }
  )
  T_ <- max(sequence_lengths)
  obs <- lapply(seq_len(n_channels), function(i) {
    suppressWarnings(suppressMessages(
      seqdef(matrix(symbol_names[[i]][1], n_sequences, T_),
             alphabet = symbol_names[[i]]
      )))
  })
  names(obs) <- paste0("Channel ", seq_len(n_channels))
  model <- build_nhmm(
    obs, n_states, initial_formula, transition_formula, emission_formula, 
    data, time, id
  )
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
  model$etas <- create_initial_values(coefs, model, init_sd)
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  if (n_channels == 1) {
    model$gammas$B <- eta_to_gamma_cube(model$etas$B)
  } else {
    model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
  }
  if (n_channels == 1L) {
    out <- simulate_nhmm_singlechannel(
      model$etas$pi, model$X_pi, 
      model$etas$A, model$X_A, 
      model$etas$B, model$X_B
    )
  } else {
    out <- simulate_nhmm_multichannel(
      model$etas$pi, model$X_pi, 
      model$etas$A, model$X_A, 
      model$etas$B, model$X_B, 
      model$n_symbols
    )
  }
  for (i in seq_len(model$n_sequences)) {
    Ti <- sequence_lengths[i]
    if (Ti < T_) {
      out$states[(Ti + 1):T_, i] <- NA
      out$observations[, (Ti + 1):T_, i] <- NA
    }
  }
  state_names <- model$state_names
  symbol_names <- model$symbol_names
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
  
  if (n_channels == 1) {
    dim(out$observations) <- dim(out$observations)[2:3]
    out$observations[] <- symbol_names[c(out$observations) + 1]
    model$observations <- suppressWarnings(suppressMessages(
      seqdef(t(out$observations), alphabet = symbol_names, cnames = seq_len(T_))
    ))
  } else {
    model$observations <- lapply(seq_len(n_channels), function(i) {
      out$observations[i, , ] <- symbol_names[[i]][c(out$observations[i, , ]) + 1]
      suppressWarnings(suppressMessages(
        seqdef(t(out$observations[i, , ]), alphabet = symbol_names[[i]], 
               cnames = seq_len(T_))
      ))
    })
    names(model$observations) <- model$channel_names
  }
  list(model = model, states = states)
}
