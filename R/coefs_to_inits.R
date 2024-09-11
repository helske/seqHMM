coefs_to_inits_mnhmm_initial <- function(x) {
  variables <- unique(x$variable)
  clusters <- unique(x$cluster)
  states <- unique(x$state)
  
  K <- length(variables)
  D <- length(clusters)
  S <- length(states)
  
  z <- array(0, c(D, S, K))
  for (i in 1:nrow(x)) {
    var_index <- match(x$variable[i], states)
    cluster_index <- match(x$cluster[i], clusters)
    state_index <- match(x$state[i], states)
    z[cluster_index, state_index, var_index] <- x$value[i]
  }
  z
}
coefs_to_inits_mnhmm_transition <- function(x) {
  variables <- unique(x$variable)
  clusters <- unique(x$cluster)
  state_from <- unique(x$state_from)
  state_to <- unique(x$state_to)
  
  K <- length(variables)
  D <- length(clusters)
  S <- length(state_from)
  
  z <- array(0, c(D, S, S - 1, K))
  for (i in 1:nrow(x)) {
    var_index <- match(x$variable[i], states)
    cluster_index <- match(x$cluster[i], clusters)
    state_from_index <- match(x$state_from[i], state_from)
    state_to_index <- match(x$state_to[i], state_to)
    z[cluster_index, state_from_index, state_to_index, var_index] <- x$value[i]
  }
  z
}