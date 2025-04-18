#' Simulate Mixture Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' mixture non-homogeneous hidden Markov model.
#'
#' @param n_sequences The number of sequences to simulate.
#' @param sequence_lengths The lengths of the simulated sequences. 
#' Either a scalar or vector of length `n_sequences`.  
#' @param n_symbols A scalar or vector of length `n_channels` giving the number 
#' of observed symbols per channel.
#' @param n_states The number of hidden states.
#' @param n_clusters The number of clusters/mixtures.
#' @param coefs If `NULL` (default), model parameters `eta_pi`, `eta_A`, 
#' `eta_B`, and `eta_omega` are generated randomly. You can also specify the coefficients by 
#' providing a named list with elements `eta_pi`, `eta_A`, `eta_B`, `eta_omega`. 
#' Alternatively, values `initial_probs`, `transition_probs`,
#' `emission_probs`, and `cluster_probs` can be used, which define the values of 
#' the \eqn{\eta} coefficients corresponding to the intercept terms, while rest 
#' of the values are generated randomly.
#' @param init_sd Standard deviation of the normal distribution used to generate
#' random coefficient values. Default is `2`.
#' @inheritParams estimate_mnhmm
#'
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_mnhmm <- function(
    n_sequences, sequence_lengths, n_symbols, n_states, n_clusters,
    initial_formula = ~1, transition_formula = ~1, emission_formula = ~1, 
    cluster_formula = ~1, data, id, time, responses, coefs = NULL, 
    init_sd = 2) {
  
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
    !missing(n_clusters) && checkmate::test_int(x = n_clusters, lower = 2L), 
    "Argument {.arg n_clusters} must be a single positive integer larger than 1."
  )
  stopifnot_(
    !missing(data),
    "{.arg data}  is missing, use {.fn simulate_mhmm} instead."
  )
  stopifnot_(
    checkmate::test_character(responses),
    "Argument {.arg responses} must be a scalar or a vector of type character."
  )
  stopifnot_(
    all(idx <- !(responses %in% names(data))),
    "Variable{?s} with name{?s} {responses[!idx]} found in {.arg data}.
    Use other name{?s} for the response variable{?s}."
  )
  stopifnot_(
    length(n_symbols) == length(responses),
    "`length(symbols) should match `length(responses)`."
  )
  sequence_lengths <- rep(sequence_lengths, length.out = n_sequences)
  n_channels <- length(n_symbols)
  T_ <- max(sequence_lengths)
  
  for (i in seq_len(n_channels)) {
    data[[responses[i]]] <- factor(1, levels = seq_len(n_symbols[i]))
  }
  model <- build_mnhmm(
    responses, n_states, n_clusters, initial_formula, transition_formula, 
    emission_formula, cluster_formula, data, id, time, scale = TRUE
  )
  if (identical(coefs, "random")) {
    coefs <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL,
      cluster_probs = NULL)
  } else {
    if (is.null(coefs$initial_probs)) coefs$initial_probs <- NULL
    if (is.null(coefs$transition_probs)) coefs$transition_probs <- NULL
    if (is.null(coefs$emission_probs)) coefs$emission_probs <- NULL
    if (is.null(coefs$cluster_probs)) coefs$cluster_probs <- NULL
  }
  model$etas <- stats::setNames(
    create_initial_values(coefs, model, init_sd), c("pi", "A", "B", "omega")
  )
  model$gammas$pi <- c(eta_to_gamma_mat_field(
    model$etas$pi
  ))
  model$gammas$A <- c(eta_to_gamma_cube_field(
    model$etas$A
  ))
  if (model$n_channels == 1L) {
    model$gammas$B <- c(eta_to_gamma_cube_field(
      model$etas$B
    ))
  } else {
    l <- lengths(model$etas$B)
    gamma_B <- c(eta_to_gamma_cube_field(unlist(model$etas$B, recursive = FALSE)))
    model$gammas$B <- split(gamma_B, rep(seq_along(l), l))
  }
  model$gammas$omega <- eta_to_gamma_mat(
    model$etas$omega
  )
  if (n_channels == 1L) {
    out <- simulate_mnhmm_singlechannel(
      model$etas$pi, model$X_pi, 
      model$etas$A, model$X_A, 
      model$etas$B, model$X_B,
      model$etas$omega, model$X_omega
    )
  } else {
    eta_B <- unlist(model$etas$B, recursive = FALSE)
    out <- simulate_mnhmm_multichannel(
      model$etas$pi, model$X_pi, 
      model$etas$A, model$X_A, 
      eta_B, model$X_B,
      model$etas$omega, model$X_omega,
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
  for (i in seq_len(n_channels)) {
    model$data[[responses[i]]] <- factor(c(out$observations[i, , ]) + 1L)
  }
  state_names <- paste0(
    rep(model$cluster_names, each = model$n_states), 
    ": ", unlist(model$state_names)
  )
  states <- cbind(
    model$data[, list(id, time), env = list(id = id, time = time)], 
    state = state_names[c(out$states) + 1]
  )
  list(model = model, states = states)
}
