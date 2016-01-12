#' Simulate Parameters of Hidden Markov Models
#' 
#' These are helper functions for quick construction of initial values for various 
#' model building functions. 
#' Mostly useful for global optimization algorithms which do not depend on initial values.
#' 
#' 
#' @export
#' @param n_states Number of states in each cluster.
#' @param n_clusters Number of clusters.
#' @param left_right Constrain the transition probabilities to upper triangular. 
#' Default is \code{FALSE}.
#' @param diag_c In case of left right model, a constant value to be added to diagonal of transition matrices before scaling.
#' @param n_symbols Number of distinct symbols in each channel.
#' @rdname simulate_pars
#' @seealso \code{\link{build_hmm}}, \code{\link{build_mhmm}},
#'   \code{\link{build_mm}}, \code{\link{build_mmm}}, and  \code{\link{build_lcm}}
#'   for constructing different types of models.
simulate_initial_probs <- function(n_states, n_clusters = 1){
  
  n_states <- rep(n_states, length = n_clusters)
  
  if(n_clusters == 1){
    x <- runif(n_states)
    x/sum(x)
  } else {
    probs <- vector("list", n_clusters)
    for(i in 1:n_clusters){
      x <- runif(n_states[i])
      probs[[i]] <- x/sum(x)
    }
    probs
  }
}
#' @export
#' @rdname simulate_pars
simulate_transition_probs <- function(n_states, n_clusters = 1, left_right = FALSE, diag_c = 0){
  
  n_states <- rep(n_states, length = n_clusters)
  if (n_clusters == 1) {
    x <- matrix(runif(n_states ^ 2), n_states, n_states) + if(left_right) diag(diag_c, n_states) else 0
    if(left_right) x[lower.tri(x)] <- 0
    probs <- x/rowSums(x)
  } else {
    probs <- vector("list", n_clusters)
    for(i in 1:n_clusters){
      x <- matrix(runif(n_states[i] ^ 2), n_states[i], n_states[i]) + if(left_right) diag(diag_c, n_states[i]) else 0
      if(left_right) x[lower.tri(x)] <- 0
      probs[[i]] <- x/rowSums(x)
    }
  }
  probs
}
#' @export
#' @rdname simulate_pars
simulate_emission_probs <- function(n_states, n_symbols, n_clusters = 1){
  n_channels <- length(n_symbols)
  emiss <- vector("list", n_clusters)
  n_states <- rep(n_states, length = n_clusters)
  if (n_channels > 1) {
    for (i in 1:n_clusters) {
      emiss[[i]] <- vector("list", n_channels)
      for (j in 1:n_channels) {
        emiss[[i]][[j]] <- matrix(runif(n_states[i] * n_symbols[j]), n_states[i], n_symbols[j])
        emiss[[i]][[j]] <- emiss[[i]][[j]]/rowSums(emiss[[i]][[j]])
      }
    }
  } else {
    for (i in 1:n_clusters) {
      emiss[[i]] <- matrix(runif(n_states[i] * n_symbols), n_states[i], n_symbols)
      emiss[[i]] <- emiss[[i]] / rowSums(emiss[[i]])
    }
  }
  if (n_clusters == 1) {
    emiss[[1]]
  } else {
    emiss
  }
}