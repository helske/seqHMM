simulate_initial_probs <- function(n_states, n_clusters = 1){
  
  drop(replicate(n_clusters, {
    x <- runif(n_states)
    x/sum(x)
  }, simplify = (n_clusters == 1)))
  
}

simulate_transition_probs <- function(n_states, n_clusters = 1, left_right = FALSE, dI = rep(0, n_states)){
  dI <- rep(dI, length = n_states)
  if (n_clusters == 1) {
    x <- matrix(runif(n_states ^ 2), n_states, n_states) + if(left_right) diag(dI) else 0
    if(left_right) x[lower.tri(x)] <- 0
    x/rowSums(x)
  } else {
    replicate(n_clusters, {
      x <- matrix(runif(n_states ^ 2), n_states, n_states) + if(left_right) diag(dI) else 0
      if(left_right) x[lower.tri(x)] <- 0
      x/rowSums(x)
    }, simplify = FALSE)
  }
}

simulate_emission_probs <- function(n_states, n_symbols, n_clusters = 1){
  n_channels <- length(n_symbols)
  emiss <- vector("list", n_clusters)
  if(n_channels > 1){
    for(i in 1:n_clusters){
      emiss[[i]] <- vector("list", n_channels)
      for(j in 1:n_channels){
        emiss[[i]][[j]] <- matrix(runif(n_states * n_symbols[j]), n_states, n_symbols[j])
        emiss[[i]][[j]] <- emiss[[i]][[j]]/rowSums(emiss[[i]][[j]])
      }
    }
  } else {
    for(i in 1:n_clusters){
      emiss[[i]] <- matrix(runif(n_states * n_symbols), n_states, n_symbols)
      emiss[[i]] <- emiss[[i]]/rowSums(emiss[[i]])
    }
  }
  if (n_clusters == 1) {
    emiss[[1]]
  } else {
    emiss
  }
}