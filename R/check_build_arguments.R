#' Create observations for the model objects
#' 
#'@noRd
#'
check_observations <- function(observations, channel_names = NULL) {
  multichannel <- is_multichannel(observations)
  # Single channel but observations is a list
  if (is.list(observations) && 
      !inherits(observations, "stslist") && 
      length(observations) == 1) {
    observations <- observations[[1]]
    multichannel <- FALSE
  }
  n_channels <- ifelse(multichannel, length(observations), 1L)
  if (multichannel) {
    n_sequences <- nrow(observations[[1]])
    length_of_sequences <- ncol(observations[[1]])
    
    symbol_names <- lapply(observations, alphabet)
    n_symbols <- lengths(symbol_names)
    
    if (length(unique(sapply(observations, nrow))) > 1L) {
      stop("The number of subjects (rows) is not the same in all channels.")
    }
    if (length(unique(sapply(observations, ncol))) > 1L) {
      stop(paste0("The length of the sequences (number of columns) is not the ",
                  "same in all channels."))
    }
    if (is.null(channel_names)) {
      if (is.null(channel_names <- names(observations))) {
        channel_names <- paste("Channel", seq_len(n_channels))
      }
    } else if (length(channel_names) != n_channels) {
      warning(paste0("The length of argument 'channel_names' does not match the ",
                     "number of channels. Names were not used."))
      channel_names <- paste("Channel", seq_len(n_channels))
    }
  } else {
    n_sequences <- nrow(observations)
    length_of_sequences <- ncol(observations)
    
    symbol_names <- alphabet(observations)
    n_symbols <- length(symbol_names)
    
    if (is.null(channel_names)) {
      channel_names <- "Observations"
    }
  }
  if (n_channels > 1L) {
    nobs <- sum(sapply(observations, function(x) {
      sum(!(x == attr(observations[[1]], "nr") |
              x == attr(observations[[1]], "void") |
              is.na(x)))
    })) / n_channels
  } else {
    nobs <- sum(!(observations == attr(observations, "nr") |
                    observations == attr(observations, "void") |
                    is.na(observations)))
  }
  
  attr(observations, "multichannel") <- multichannel
  attr(observations, "n_channels") <- ifelse(
    multichannel, length(observations), 1L)
  attr(observations, "n_sequences") <- n_sequences
  attr(observations, "length_of_sequences") <- length_of_sequences
  attr(observations, "symbol_names") <- symbol_names
  attr(observations, "n_symbols") <- n_symbols
  attr(observations, "channel_names") <- channel_names
  attr(observations, "nobs") <- nobs
  observations
}

check_transition_probs <- function(transition_probs, state_names = NULL) {
  if (!is.matrix(transition_probs)) {
    stop(paste("Object provided for 'transition_probs' is not a matrix."))
  }
  if (dim(transition_probs)[1] != dim(transition_probs)[2]) {
    stop("Object 'transition_probs' must be a square matrix.")
  }
  n_states <- nrow(transition_probs)
  if (is.null(state_names)) {
    if (is.null(state_names <- rownames(transition_probs))) {
      state_names <- paste("State", seq_len(n_states))
    }
  } else {
    if (length(state_names) != n_states) {
      stop("Length of 'state_names' is not equal to the number of hidden states.")
    }
  }
  if (!isTRUE(all.equal(
    rowSums(transition_probs), 
    rep(1, dim(transition_probs)[1]), check.attributes = FALSE))
  ) {
    stop("Transition probabilities in 'transition_probs' do not sum to one.")
  }
  dimnames(transition_probs) <- list(from = state_names, to = state_names)
  transition_probs
}

check_initial_probs <- function(initial_probs, n_states, state_names = NULL) {
  if (!is.vector(initial_probs)) {
    stop(paste("Object provided for 'initial_probs' is not a vector."))
  }
  if (length(initial_probs) != n_states) {
    stop(paste("Length of 'initial_probs' is not equal to the number of states."))
  }
  if (!isTRUE(all.equal(sum(initial_probs), 1, check.attributes = FALSE))) {
    stop("Initial state probabilities in 'initial_probs' do not sum to one.")
  }
  names(initial_probs) <- state_names
  initial_probs
}

check_emission_probs <- function(
    emission_probs, n_states, n_channels, n_symbols, state_names, symbol_names,
    channel_names = NULL) {
  
  if (is.list(emission_probs) && length(emission_probs) == 1L) {
    emission_probs <- emission_probs[[1]]
  }
  if (is.list(emission_probs)) {
    if (length(emission_probs) != n_channels) {
      stop(paste0("Number of channels defined by 'emission_probs' differs ",
                  "from one defined by observations."))
    }
    for (j in seq_len(n_channels)) {
      if (!is.matrix(emission_probs[[j]])) {
        stop(paste0("Object provided in 'emission_probs' for channel ", j, 
                    " is not a matrix."))
      }
      if (!isTRUE(all.equal(
        rowSums(emission_probs[[j]]), 
        rep(1, n_states), 
        check.attributes = FALSE))) {
        stop("Emission probabilities in 'emission_probs' for channel ", j, 
             "do not sum to one.")
      }
      if (nrow(emission_probs[[j]]) != n_states) {
        stop(paste0("Number of rows in 'emission_probs' for channel ", j, 
                    " is not equal to the number of states."))
      }
      if (ncol(emission_probs[[j]]) != n_symbols[j]) {
        stop(paste0("Number of columns in 'emission_probs' for channel ", j, 
                    " is not equal to the number of symbols."))
      }
      
    }
    if (is.null(channel_names)) {
      if (is.null(channel_names <- names(observations))) {
        channel_names <- paste("Channel", 1:n_channels)
      }
    } else if (length(channel_names) != n_channels) {
      warning(paste0("The length of argument 'channel_names' does not match ", 
                     "the number of channels. Names were not used."))
      channel_names <- paste("Channel", 1:n_channels)
    }
    for (i in seq_len(n_channels)) {
      dimnames(emission_probs[[i]]) <- 
        list(state_names = state_names, symbol_names = symbol_names[[i]])
    }
    names(emission_probs) <- channel_names
  } else {
    if (!is.matrix(emission_probs)) {
      stop(paste("Object provided for 'emission_probs' is not a matrix."))
    }
    if (n_states != nrow(emission_probs)) {
      stop("Number of rows in 'emission_probs' is not equal to the number of states.")
    }
    if (n_symbols != ncol(emission_probs)) {
      stop("Number of columns in 'emission_probs' is not equal to the number of symbols.")
    }
    if (!isTRUE(all.equal(rowSums(emission_probs), rep(1, n_states), 
                          check.attributes = FALSE))) {
      stop("Emission probabilities in 'emission_probs' do not sum to one.")
    }
    dimnames(emission_probs) <- 
      list(state_names = state_names, symbol_names = symbol_names)
  }
  emission_probs
}