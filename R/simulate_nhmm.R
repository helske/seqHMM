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
  symbol_names <- lapply(seq_len(n_channels), function(i) seq_len(n_symbols[i]))
  obs <- lapply(seq_len(n_channels), function(i) {
    suppressWarnings(suppressMessages(
      seqdef(matrix(symbol_names[[i]][1], n_sequences, max(sequence_lengths)),
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
  K_i <- nrow(model$X_initial)
  K_s <- nrow(model$X_transition)
  K_o <- nrow(model$X_emission)
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  if (n_channels == 1L) {
    model$gammas$B <- eta_to_gamma_cube(model$etas$B)
  } else {
    model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
  }
  probs <- get_probs(model)
  states <- array(NA_character_, c(max(sequence_lengths), n_sequences))
  obs <- array(NA_character_, c(max(sequence_lengths), n_channels, n_sequences))
  ids <- unique(data[[id]])
  times <- sort(unique(data[[time]]))
  state_names <- model$state_names
  for (i in seq_len(n_sequences)) {
    p_init <- probs$initial[
      probs$initial[[time]] == time[1] & probs$initial[[id]] == ids[i],
      "probability"
    ]
    states[1, i] <- sample(state_names, 1, prob = p_init)
    for (k in seq_len(n_channels)) {
      p_emission <- probs$emission[
        probs$emission[[time]] == time[1] & probs$emission[[id]] == ids[i] &
          probs$emission$state == states[1, i] & probs$emission$channel == k,
        "probability"
      ]
      obs[1, k, i] <- sample(symbol_names[[k]], 1, prob = p_emission)
    }
  }
  
  for (i in seq_len(n_sequences)) {
    for (t in 2:sequence_lengths[i]) {
      p_transition <- probs$transition[
        probs$transition[[time]] == times[t] & probs$transition[[id]] == ids[i] &
          probs$transition$state_from == states[t - 1, i], "probability"
      ]
      states[t, i] <- sample(state_names, 1, prob = p_transition)
      for (k in seq_len(n_channels)) {
        p_emission <- probs$emission[
          probs$emission[[time]] == time[t] & probs$emission[[id]] == ids[i] &
            probs$emission$state == states[t, i] & probs$emission$channel == k,
          "probability"
        ]
        obs[t, k, i] <- sample(symbol_names[[k]], 1, prob = p_emission)
      }
    }
  }
  states <- suppressWarnings(suppressMessages(
    seqdef(
      matrix(
        t(states),
        n_sequences, max(sequence_lengths)
      ), 
      alphabet = state_names
    )
  ))
  obs <- lapply(seq_len(n_channels), function(i) {
    suppressWarnings(suppressMessages(
      seqdef(t(obs[, i, ]), alphabet = symbol_names[[i]])
    ))
  })
  names(obs) <- model$channel_names
  if (n_channels == 1) obs <- obs[[1]]
  model$observations <- obs
  list(model = model, states = states)
}
