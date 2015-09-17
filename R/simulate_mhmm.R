#' Simulate Mixture Hidden Markov Models
#' 
#' Simulate sequences of observed and hidden states given parameters of a mixture hidden Markov model.
#'
#' @param n_sequences Number of simulations.
#' @param initial_probs A list containing vectors of initial state probabilities 
#' for submodels of each cluster.
#' @param transition_matrix A list of matrices of transition probabilities for submodels of each cluster.
#' @param emission_matrix A list which contains matrices of emission probabilities or a list of such 
#' objects (one for each channel) for submodels of each cluster. Note that the matrices must have 
#' dimensions m x s where m is the number of hidden states and s is the number of unique symbols 
#' (observed states) in the data.
#' @param sequence_length Length for simulated sequences.
#' @param formula Covariates as an object of class \code{\link{formula}}, left side omitted.
#' @param data An optional data frame, list or environment containing the variables in the model. 
#' If not found in data, the variables are taken from \code{environment(formula)}.
#' @param An optional k x l matrix of regression coefficients for time-constant covariates 
#' for mixture probabilities, where l is the number of clusters and k is the number of 
#' covariates. A logit-link is used for mixture probabilities. The first column is set to zero.
#' 
#' @return A list of state sequence objects of class \code{stslist}.
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_mhmm}} for building 
#' and fitting mixture hidden Markov models; \code{\link{ssplot}} for plotting 
#' multiple sequence data sets; \code{\link{seqdef}} for more
#' information on state sequence objects; and \code{\link{simulate_hmm}}
#' for simulating hidden Markov models.
#' @export
simulate_mhmm <- function(n_sequences, initial_probs, transition_matrix, 
  emission_matrix, sequence_length, formula, data, coef){
  
  if (is.list(transition_matrix)){
    n_clusters<-length(transition_matrix)
  } else {
    stop("transition_matrix is not a list.")
  }
  if (length(emission_matrix)!=n_clusters || length(initial_probs)!=n_clusters) {
    stop("Unequal list lengths of transition_matrix, emission_matrix and initial_probs.")
  }
  if (is.null(cluster_names <- names(transition_matrix))) {
    cluster_names <- paste("Cluster", 1:n_clusters)
  }
  
  if (is.list(emission_matrix[[1]])) {
    n_channels <- length(emission_matrix[[1]])
  } else {
    n_channels <- 1
    for(i in 1:n_clusters)
      emission_matrix[[i]] <- list(emission_matrix[[i]])
  }
  if (is.null(channel_names <- names(emission_matrix[[1]]))) {
    channel_names <- 1:n_channels
  }
  if (n_sequences < 2) {
    stop("Number of simulations (n_sequences) must be at least 2 for a mixture model.")
  }
  
  if (missing(formula)) {
    formula <- stats::formula(rep(1, n_sequences) ~ 1)
  }
  if (missing(data)) {
    data <- environment(formula)
  }
  if (inherits(formula, "formula")) {
    X <- model.matrix(formula, data)
    if (nrow(X) != n_sequences) {
      if (length(all.vars(formula)) > 0 && 
          sum(!complete.cases(data[all.vars(formula)])) > 0) {
        stop("Missing cases are not allowed in covariates. Use e.g. the complete.cases function to detect them, then fix, impute, or remove.") 
      } else {
        stop("Number of subjects in data for covariates does not match the number of subjects in the sequence data.")
      }
    }
    n_covariates <- ncol(X)
  } else {
    stop("Object given for argument formula is not of class formula.")
  }
  if (missing(coef)) {
    coef <- matrix(0, n_covariates, n_clusters)
  } else {
    if (ncol(coef) != n_clusters | nrow(coef) != n_covariates) {
      stop("Wrong dimensions of coef.")
    }
    coef[, 1] <- 0
  }       
  
  pr <- exp(X %*% coef)
  pr <- pr / rowSums(pr)
  
  
  n_symbols <- sapply(emission_matrix[[1]], ncol)
  if (is.null(colnames(emission_matrix[[1]][[1]]))) {
    symbol_names <- lapply(1:n_channels, function(i) 1:n_symbols[i])
  } else symbol_names <- lapply(1:n_channels, function(i) colnames(emission_matrix[[1]][[i]]))
  
  obs <- lapply(1:n_channels, function(i) suppressWarnings(suppressMessages(seqdef(matrix(NA, n_sequences, sequence_length), 
    alphabet = symbol_names[[i]]))))
  
  names(obs) <- channel_names
  
  n_states <- sapply(transition_matrix, nrow)
  if (is.null(rownames(transition_matrix[[1]]))) {
    state_names <- lapply(1:n_clusters, function(i) 1:n_states[i])
  } else state_names <- lapply(1:n_clusters, function(i) rownames(transition_matrix[[i]]))
  v_state_names <- unlist(state_names)
  if (length(unique(v_state_names)) != length(v_state_names)){
    for (i in 1:n_clusters) {
      colnames(transition_matrix[[i]]) <- rownames(transition_matrix[[i]]) <- 
        paste(cluster_names[i], state_names[[i]], sep=":")
    }
  }
  for (i in 1:n_clusters) {
    for (j in 1:n_channels) {
      rownames(emission_matrix[[i]][[j]]) <- 
        paste(cluster_names[i], state_names[[i]], sep=":")
    }
  }
  
  v_state_names <- paste(rep(cluster_names, n_states), v_state_names, sep=":")
  
  states <- suppressWarnings(suppressMessages(seqdef(matrix(NA, 
    n_sequences, sequence_length), alphabet = v_state_names)))
  
  clusters <- numeric(n_sequences)
  for (i in 1:n_sequences) {
    clusters[i] <- sample(cluster_names, size = 1, prob = pr[i, ])
  }
  for (i in 1:n_clusters) {
    if(sum(clusters == cluster_names[i]) > 0) {
      sim <- simulate_hmm(n_sequences = sum(clusters == cluster_names[i]), initial_probs[[i]],
        transition_matrix[[i]], emission_matrix[[i]], sequence_length)
      for (k in 1:n_channels) {
        obs[[k]][clusters == cluster_names[i], ] <- sim$observations[[k]]
      }
      states[clusters == cluster_names[i], ] <- sim$states
    }
  }
  

  p <- 0
  for (i in 1:n_channels) {
    attr(obs[[i]], "cpal") <- seqHMM::colorpalette[[
      length(unlist(symbol_names))]][(p + 1):(p + n_symbols[[i]])]
    p <- p + n_symbols[[i]]
  }
  
  if (length(unlist(symbol_names)) != length(alphabet(states))) {
    attr(states, "cpal") <- seqHMM::colorpalette[[length(alphabet(states))]]
  } else {
    attr(states, "cpal") <- seqHMM::colorpalette[[length(alphabet(states)) + 1]][1:length(alphabet(states))]
  }
  
  
  list(observations = obs, states = states)
}