#' Simulate Parameters of Hidden Markov Models
#'
#' These are helper functions for quick construction of initial values for various
#' model building functions.
#' Mostly useful for global optimization algorithms which do not depend on initial values.
#'
#' @export
#' @param n_states Number of states in each cluster.
#' @param n_clusters Number of clusters.
#' @param left_right Constrain the transition probabilities to upper triangular.
#' Default is `FALSE`.
#' @param diag_c A constant value to be added to diagonal of transition matrices before scaling.
#' @param n_symbols Number of distinct symbols in each channel.
#' @param alpha A scalar, or a vector of length S (number of states) or M 
#' (number of symbols) defining the parameters of the Dirichlet distribution 
#' used to simulate the probabilities.
#' @rdname simulate_pars
simulate_initial_probs <- function(n_states, n_clusters = 1, alpha = 1) {
  n_states <- rep(n_states, length = n_clusters)
  
  if (n_clusters == 1) {
    x <- stats::rgamma(n_states, alpha)
    x / sum(x)
  } else {
    probs <- vector("list", n_clusters)
    for (i in 1:n_clusters) {
      x <- stats::rgamma(n_states[i], alpha)
      probs[[i]] <- x / sum(x)
    }
    probs
  }
}
#' @export
#' @rdname simulate_pars
simulate_transition_probs <- function(n_states, n_clusters = 1, 
                                      left_right = FALSE, diag_c = 0, alpha = 1) {
  n_states <- rep(n_states, length = n_clusters)
  if (n_clusters == 1) {
    x <- matrix(stats::rgamma(n_states^2, alpha), n_states, n_states, TRUE) + 
      diag(diag_c, n_states)
    if (left_right) x[lower.tri(x)] <- 0
    probs <- x / rowSums(x)
  } else {
    probs <- vector("list", n_clusters)
    for (i in seq_len(n_clusters)) {
      x <- matrix(
        stats::rgamma(n_states[i]^2, alpha), n_states[i], n_states[i], TRUE
      ) + diag(diag_c, n_states[i])
      if (left_right) x[lower.tri(x)] <- 0
      probs[[i]] <- x / rowSums(x)
    }
  }
  probs
}
#' @export
#' @rdname simulate_pars
simulate_emission_probs <- function(n_states, n_symbols, n_clusters = 1, 
                                    alpha = 1) {
  n_channels <- length(n_symbols)
  emiss <- vector("list", n_clusters)
  n_states <- rep(n_states, length = n_clusters)
  if (n_channels > 1) {
    for (i in seq_len(n_clusters)) {
      emiss[[i]] <- vector("list", n_channels)
      for (j in seq_len(n_channels)) {
        emiss[[i]][[j]] <- matrix(
          stats::rgamma(n_states[i] * n_symbols[j], alpha), n_states[i], n_symbols[j],
          TRUE)
        emiss[[i]][[j]] <- emiss[[i]][[j]] / rowSums(emiss[[i]][[j]])
      }
    }
  } else {
    for (i in seq_len(n_clusters)) {
      emiss[[i]] <- matrix(
        stats::rgamma(n_states[i] * n_symbols, alpha), n_states[i], n_symbols, 
        TRUE
      )
      emiss[[i]] <- emiss[[i]] / rowSums(emiss[[i]])
    }
  }
  if (n_clusters == 1) {
    emiss[[1]]
  } else {
    emiss
  }
}
