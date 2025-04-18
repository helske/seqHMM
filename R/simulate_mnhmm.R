#' Simulate Mixture Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' mixture non-homogeneous hidden Markov model.
#'
#' @param n_clusters The number of clusters/mixtures.
#' @param cluster_formula of class [formula()] for the mixture probabilities.
#' @param coefs Same as argument `inits` in [estimate_mnhmm()]. If `NULL` 
#' (default), the model parameters are generated randomly.
#' @inheritParams simulate_nhmm
#'
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_mnhmm <- function(
    n_states, n_clusters, emission_formula, initial_formula = ~1, 
    transition_formula = ~1, cluster_formula = ~ 1, data, id, time, 
    coefs = NULL, init_sd = 2) {
  
  stopifnot_(
    !missing(data),
    "{.arg data} is missing, use {.fn simulate_mhmm} instead."
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
  cluster_names <- paste("Cluster", seq_len(n_clusters))
  if (is.null(coefs)) {
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
  model <- build_mnhmm(
    n_states, n_clusters, emission_formula, initial_formula, transition_formula, 
    cluster_formula, data, id, time, coefs = coefs
  )
  out <- Rcpp_simulate_mnhmm(
    create_obsArray(model), model$sequence_lengths, model$n_symbols, 
    model$X_pi, model$X_A, model$X_B, model$X_omega,
    io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
    iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
    model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, model$gammas$gamma_omega
  )
  T_ <- model$length_of_sequences
  for (i in seq_len(model$n_sequences)) {
    Ti <- model$sequence_lengths[i]
    if (Ti < T_) {
      out$states[(Ti + 1):T_, i] <- NA
      out$observations[, (Ti + 1):T_, i] <- NA
    }
  }
  for (i in seq_len(model$n_channels)) {
    model$data[[model$responses[i]]][] <- 
      model$symbol_names[[i]][c(out$observations[i, , ] + 1)]
  }
  state_names <- paste0(
    rep(model$cluster_names, each = model$n_states), ": ",
    unlist(model$state_names)
  )
  states <- cbind(
    model$data[, list(id, time), env = list(id = id, time = time)], 
    state = state_names[c(out$states) + 1]
  )
  list(model = model, states = states)
}
