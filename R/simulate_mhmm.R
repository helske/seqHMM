#' Simulate Mixture Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a mixture
#' hidden Markov model.
#'
#' @param n_sequences The number of sequences to simulate.
#' @param initial_probs A list containing vectors of initial state probabilities
#' for the submodel of each cluster.
#' @param transition_probs A list of matrices of transition probabilities
#'   for the submodel of each cluster.
#' @param emission_probs A list which contains matrices of emission
#'   probabilities or a list of such objects (one for each channel) for
#'   the submodel of each cluster. Note that the matrices must have
#'   dimensions \eqn{s x m} where \eqn{s} is the number of hidden states
#'   and \eqn{m} is the number of unique symbols (observed states) in the
#'   data.
#' @param sequence_length The length of the simulated sequences.
#' @param formula Covariates as an object of class [formula()],
#'   left side omitted.
#' @param data An optional data frame, a list or an environment containing
#'   the variables in the model. If not found in data, the variables are
#'   taken from `environment(formula)`.
#' @param coefficients An optional \eqn{k x l} matrix of regression
#'   coefficients for time-constant covariates for mixture probabilities,
#'   where \eqn{l} is the number of clusters and \eqn{k} is the number of
#'   covariates. A logit-link is used for mixture probabilities. The first
#'   column is set to zero.
#'
#' @return A list of state sequence objects of class `stslist`.
#' @seealso [build_mhmm()] and [fit_model()] for building
#' and fitting mixture hidden Markov models.
#' @export
#' @examples
#' emission_probs_1 <- matrix(c(0.75, 0.05, 0.25, 0.95), 2, 2)
#' emission_probs_2 <- matrix(c(0.1, 0.8, 0.9, 0.2), 2, 2)
#' colnames(emission_probs_1) <- colnames(emission_probs_2) <-
#'   c("heads", "tails")
#'
#' transition_probs_1 <- matrix(c(9, 0.1, 1, 9.9) / 10, 2, 2)
#' transition_probs_2 <- matrix(c(35, 1, 1, 35) / 36, 2, 2)
#' rownames(emission_probs_1) <- rownames(transition_probs_1) <-
#'   colnames(transition_probs_1) <- c("coin 1", "coin 2")
#' rownames(emission_probs_2) <- rownames(transition_probs_2) <-
#'   colnames(transition_probs_2) <- c("coin 3", "coin 4")
#'
#' initial_probs_1 <- c(1, 0)
#' initial_probs_2 <- c(1, 0)
#'
#' n <- 30
#' set.seed(123)
#' covariate_1 <- runif(n)
#' covariate_2 <- sample(c("A", "B"),
#'   size = n, replace = TRUE,
#'   prob = c(0.3, 0.7)
#' )
#' dataf <- data.frame(covariate_1, covariate_2)
#'
#' coefs <- cbind(cluster_1 = c(0, 0, 0), cluster_2 = c(-1.5, 3, -0.7))
#' rownames(coefs) <- c("(Intercept)", "covariate_1", "covariate_2B")
#'
#' sim <- simulate_mhmm(
#'   n = n, initial_probs = list(initial_probs_1, initial_probs_2),
#'   transition_probs = list(transition_probs_1, transition_probs_2),
#'   emission_probs = list(emission_probs_1, emission_probs_2),
#'   sequence_length = 20, formula = ~ covariate_1 + covariate_2,
#'   data = dataf, coefficients = coefs
#' )
#'
#' stacked_sequence_plot(sim, 
#'   sort_by = "start", sort_channel = "states", type = "i"
#' )
#'
#' hmm <- build_mhmm(sim$observations,
#'   initial_probs = list(initial_probs_1, initial_probs_2),
#'   transition_probs = list(transition_probs_1, transition_probs_2),
#'   emission_probs = list(emission_probs_1, emission_probs_2),
#'   formula = ~ covariate_1 + covariate_2,
#'   data = dataf
#' )
#'
#' fit <- fit_model(hmm)
#' fit$model
#'
#' paths <- hidden_paths(fit$model, as_stslist = TRUE)
#'
#' stacked_sequence_plot(
#'   list(
#'     "estimated paths" = paths, 
#'     "true (simulated)" = sim$states
#'   ),
#'   sort_by = "start",
#'   sort_channel = "true (simulated)",
#'   type = "i"
#' )
#'
simulate_mhmm <- function(
    n_sequences, initial_probs, transition_probs, emission_probs, 
    sequence_length, formula = NULL, data = NULL, coefficients = NULL) {
  
  stopifnot_(
    is.list(transition_probs),
    "{.arg transition_probs} is not a {.cls list}."
  )
  stopifnot_(
    is.list(emission_probs),
    "{.arg emission_probs} is not a {.cls list}."
  )
  stopifnot_(
    is.list(initial_probs),
    "{.arg initial_probs} is not a {.cls list}."
  )
  n_clusters <- length(transition_probs)
  stopifnot_(
    n_clusters > 1L,
    "{.arg transition_probs} is a list of length 1 leading to an
        ordinary HMM. Please use {.fn simulate_hmm} instead."
  )
  stopifnot_(
    length(emission_probs) == n_clusters,
    "Length of a {.arg emission_probs} does not match with the length of
      {.arg transition_probs}."
  )
  stopifnot_(
    length(initial_probs) == n_clusters,
    "Length of a {.arg initial_probs} does not match with the length of
      {.arg transition_probs}."
  )
  
  if (is.list(emission_probs[[1]])) {
    n_channels <- length(emission_probs[[1]])
  } else {
    n_channels <- 1
    for (i in 1:n_clusters) {
      emission_probs[[i]] <- list(emission_probs[[i]])
    }
  }
  if (is.null(channel_names <- names(emission_probs[[1]]))) {
    channel_names <- paste("Channel", seq_len(n_channels))
  }
  if (n_sequences < 2) {
    stop_("{.arg n_sequences} must be at least 2 for a mixture model.")
  }
  n_symbols <- sapply(emission_probs[[1]], ncol)
  if (is.null(colnames(emission_probs[[1]][[1]]))) {
    symbol_names <- lapply(seq_len(n_channels), \(i) seq_len(n_symbols[i]))
  } else {
    symbol_names <- lapply(seq_len(n_channels), \(i) colnames(emission_probs[[1]][[i]]))
  }
  
  obs <- lapply(seq_len(n_channels), \(i) {
    suppressWarnings(suppressMessages(
      seqdef(matrix(symbol_names[[i]], n_sequences, sequence_length),
             alphabet = symbol_names[[i]]
      )))
  })
  names(obs) <- channel_names
  n_states <- sapply(transition_probs, nrow)
  if (is.null(rownames(transition_probs[[1]]))) {
    state_names <- lapply(seq_len(n_clusters), \(i) {
      paste("State", seq_len(n_states[i]))
    })
  } else {
    state_names <- lapply(seq_len(n_clusters), \(i) {
      rownames(transition_probs[[i]])
    }
    )
  }
  if (is.null(cluster_names <- names(transition_probs))) {
    cluster_names <- paste0("Cluster ", seq_len(n_clusters))
  }
  v_state_names <- unlist(state_names)
  if (n_unique(v_state_names) != length(v_state_names)) {
    for (i in seq_len(n_clusters)) {
      transition_probs[[i]] <- 
        .check_transition_probs(
          transition_probs[[i]], 
          paste(cluster_names[i], state_names[[i]], sep = ":"))
    }
    v_state_names <- paste(rep(cluster_names, n_states), v_state_names, sep = ":")
  } else {
    for (i in seq_len(n_clusters)) {
      transition_probs[[i]] <- 
        .check_transition_probs(transition_probs[[i]], state_names[[i]])
    }
  }
  states <- suppressWarnings(suppressMessages(
    seqdef(matrix(
      v_state_names[1],
      n_sequences, sequence_length
    ), alphabet = v_state_names)
  ))
  if (n_channels == 1) {
    for (i in 1:n_clusters) {
      emission_probs[[i]] <- emission_probs[[i]][[1]]
    }
    symbol_names <- symbol_names[[1]]
  }
  for (i in seq_len(n_clusters)) {
    initial_probs[[i]] <- 
      .check_initial_probs(initial_probs[[i]], n_states[i], state_names[[i]])
    emission_probs[[i]] <- .check_emission_probs(
      emission_probs[[i]], n_states[i], n_channels, n_symbols, 
      state_names[[i]], symbol_names
    )
  }
  
  if (is.null(formula)) {
    formula <- formula(~ 1)
    X <- stats::model.matrix(formula, data = data.frame(y = rep(1, n_sequences)))
    n_covariates <- 1L
  } else {
    stopifnot_(
      inherits(formula, "formula"), 
      "Argument {.arg formula} must be a {.cls formula} object.")
    stopifnot_(
      is.data.frame(data),
      "If {.arg formula} is provided, {.arg data} must be a {.cls data.frame} 
      object."
    )
    X <- stats::model.matrix(formula, data)
    if (nrow(X) != n_sequences) {
      if (length(all.vars(formula)) > 0 && 
          sum(!stats::complete.cases(data[all.vars(formula)])) > 0) {
        stop_(
          "Missing cases are not allowed in covariates. Use e.g. the 
                {.fn stats::complete.cases} to detect them, then fix, impute, or 
                remove."
        )
      } else {
        stop_(
          "Number of subjects in data for covariates does not match the 
            number of subjects in the sequence data."
        )
      }
    }
    n_covariates <- ncol(X)
  }
  if (is.null(coefficients)) {
    coefficients <- matrix(0, n_covariates, n_clusters)
  } else {
    if (ncol(coefficients) != n_clusters | nrow(coefficients) != n_covariates) {
      stop_(
        "Wrong dimensions of {.arg coefficients}. Should be 
        {n_clusters} x {n_covariates}"
      )
    }
    coefficients[, 1] <- 0
  }
  pr <- exp(X %*% coefficients)
  pr <- pr / rowSums(pr)
  
  clusters <- numeric(n_sequences)
  for (i in 1:n_sequences) {
    clusters[i] <- sample(cluster_names, size = 1, prob = pr[i, ])
  }
  for (i in 1:n_clusters) {
    if (sum(clusters == cluster_names[i]) > 0) {
      sim <- simulate_hmm(
        n_sequences = sum(clusters == cluster_names[i]), initial_probs[[i]],
        transition_probs[[i]], emission_probs[[i]], sequence_length
      )
      if (n_channels > 1) {
        for (k in seq_len(n_channels)) {
          obs[[k]][clusters == cluster_names[i], ] <- sim$observations[[k]]
        }
      } else {
        obs[[1]][clusters == cluster_names[i], ] <- sim$observations
      }
      states[clusters == cluster_names[i], ] <- sim$states
    }
  }
  p <- 0
  if (length(unlist(symbol_names)) <= 200) {
    for (i in seq_len(n_channels)) {
      TraMineR::cpal(obs[[i]]) <- seqHMM::colorpalette[[
        length(unlist(symbol_names))
      ]][(p + 1):(p + n_symbols[[i]])]
      p <- p + n_symbols[[i]]
    }
  } else {
    cp <- NULL
    k <- 200
    l <- 0
    while (length(unlist(symbol_names)) - l > 0) {
      cp <- c(cp, seqHMM::colorpalette[[k]])
      l <- l + k
      k <- k - 1
    }
    cp <- cp[1:length(unlist(symbol_names))]
    for (i in seq_len(n_channels)) {
      TraMineR::cpal(obs[[i]]) <- cp[(p + 1):(p + n_symbols[[i]])]
      p <- p + n_symbols[[i]]
    }
  }
  
  if (length(unlist(symbol_names)) != length(alphabet(states))) {
    if (length(alphabet(states)) <= 200) {
      TraMineR::cpal(states) <- seqHMM::colorpalette[[length(alphabet(states))]]
    } else {
      cp <- NULL
      k <- 200
      p <- 0
      while (length(alphabet(states)) - p > 0) {
        cp <- c(cp, seqHMM::colorpalette[[k]])
        p <- p + k
        k <- k - 1
      }
      TraMineR::cpal(states) <- cp[1:length(alphabet(states))]
    }
  } else {
    if (length(alphabet(states)) <= 199) {
      TraMineR::cpal(states) <- seqHMM::colorpalette[[length(alphabet(states)) + 1]][1:length(alphabet(states))]
    } else {
      cp <- NULL
      k <- 199
      p <- 0
      while (length(alphabet(states)) - p > 0) {
        cp <- c(cp, seqHMM::colorpalette[[k]])
        p <- p + k
        k <- k - 1
      }
      TraMineR::cpal(states) <- cp[1:length(alphabet(states))]
    }
  }
  if (n_channels == 1) obs <- obs[[1]]
  list(observations = obs, states = states)
}
