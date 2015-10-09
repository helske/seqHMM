#' Simulate hidden Markov models
#' 
#' Simulate sequences of observed and hidden states given parameters of a hidden Markov model.
#'
#' @param n_sequences Number of simulations.
#' @param initial_probs A vector of initial state probabilities.
#' @param transition_matrix A matrix of transition probabilities.
#' @param emission_matrix A matrix of emission probabilities or a list of such objects (one for each channel).
#' @param sequence_length Length for simulated sequences.
#'
#' @return A list of state sequence objects of class \code{stslist}.
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building 
#' and fitting hidden Markov models; \code{\link{ssplot}} for plotting 
#' multiple sequence data sets; \code{\link{seqdef}} for more
#' information on state sequence objects; and \code{\link{simulate_mhmm}}
#' for simulating mixture hidden Markov models.
#' @export
#'
#' @examples
#' # Parameters for a HMM
#' emission_matrix <- matrix(c(0.5, 0.2, 0.5, 0.8), 2, 2)
#' transition_matrix <- matrix(c(5/6, 1/6, 1/6, 5/6), 2, 2)
#' initial_probs <- c(1, 0)
#' 
#' # Setting seed for simulation
#' set.seed(1)
#' 
#' # Simulating sequences
#' sim <- simulate_hmm(
#'   n_sequences = 10, initial_probs = initial_probs, 
#'   transition_matrix = transition_matrix, 
#'   emission_matrix = emission_matrix, 
#'   sequence_length = 20)

simulate_hmm <- function(n_sequences, initial_probs, transition_matrix, emission_matrix, 
  sequence_length){
  
  if (is.list(emission_matrix)) {
    n_channels <- length(emission_matrix)
  } else {
    n_channels <- 1
    emission_matrix <- list(emission_matrix)
  }
  
  n_states <- length(initial_probs)
  if (is.null(state_names <- rownames(transition_matrix))) {
    state_names <- 1:n_states
  }
  for (i in 1:n_channels) rownames(emission_matrix[[i]]) <- rownames(transition_matrix)
  n_symbols <- sapply(emission_matrix, ncol)
  if (is.null(colnames(emission_matrix[[1]]))) {
    symbol_names <- lapply(1:n_channels, function(i) 1:n_symbols[i])
  } else symbol_names <- lapply(1:n_channels, function(i) colnames(emission_matrix[[i]]))
  
  if (is.null(channel_names <- names(emission_matrix))) {
    channel_names <- 1:n_channels
  }
  
  states <- array(NA, c(n_sequences, sequence_length))
  obs  <- array(NA, c(n_sequences, sequence_length, n_channels))
  
  for (i in 1:n_sequences) {
    states[i, 1] <- sample(state_names, 1, prob = initial_probs)
    for (k in 1:n_channels) {
      obs[i, 1, k] <- sample(symbol_names[[k]], 1, 
        prob = emission_matrix[[k]][states[i, 1], ])
    }
  }
  
  
  if (sequence_length > 1) {
    for (i in 1:n_sequences) {
      for (t in 2:sequence_length) {
        states[i, t] <- sample(state_names, 1, prob = transition_matrix[states[i, t - 1], ])
        for (k in 1:n_channels) {
          obs[i, t, k] <- sample(symbol_names[[k]], 1, 
            prob = emission_matrix[[k]][states[i, t], ])
        }
      }
    }
    
  }
  obs <- suppressMessages(lapply(1:n_channels, function(i) 
    seqdef(obs[, , i], alphabet = symbol_names[[i]])))
  
  
  states <- suppressMessages(seqdef(states, alphabet = state_names))
  
  
  p <- 0
  if (length(unlist(symbol_names)) <= 200) {
    for (i in 1:n_channels) {
      attr(obs[[i]], "cpal") <- seqHMM::colorpalette[[
        length(unlist(symbol_names))]][(p + 1):(p + n_symbols[[i]])]
      p <- 1
    }
  } else {
    cp <- NULL
    k <- 200
    l <- 0
    while(length(unlist(symbol_names)) - l > 0){
      cp <- c(cp, seqHMM::colorpalette[[k]])
      l <- l + k
      k <- k - 1
    }
    cp <- cp[1:length(unlist(symbol_names))]
    for (i in 1:n_channels) {
      attr(obs[[i]], "cpal") <- cp[(p + 1):(p + n_symbols[[i]])]
      p <- 1
    }
  }
  
  
  if (length(unlist(symbol_names)) != length(alphabet(states))) {
    if (length(alphabet(states)) <= 200) {
      attr(states, "cpal") <- seqHMM::colorpalette[[length(alphabet(states))]]
    } else {
      cp <- NULL
      k <- 200
      p <- 0
      while(length(alphabet(states)) - p > 0){
        cp <- c(cp, seqHMM::colorpalette[[k]])
        p <- p + k
        k <- k - 1
      }
      attr(states, "cpal") <- cp[1:length(alphabet(states))]
    }
  } else {
    if (length(alphabet(states)) <= 199) {
      attr(states, "cpal") <- seqHMM::colorpalette[[length(alphabet(states)) + 1]][1:length(alphabet(states))]
    } else {
      cp <- NULL
      k <- 199
      p <- 0
      while(length(alphabet(states)) - p > 0){
        cp <- c(cp, seqHMM::colorpalette[[k]])
        p <- p + k
        k <- k - 1
      }
      attr(states, "cpal") <- cp[1:length(alphabet(states))]
    }
  }
  
  if (n_channels == 1) obs <- obs[[1]]
  
  list(observations = obs, states = states)
}