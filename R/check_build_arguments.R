#' Create observations for the model objects
#' 
#'@noRd
#'
.check_observations <- function(observations, channel_names = NULL) {
  multichannel <- is_multichannel(observations)
  # Single channel but observations is a list
  if (is.list(observations) && 
      !TraMineR::is.stslist(observations) && 
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
    dim(n_symbols) <- length(n_symbols)
    if (is.null(channel_names)) {
      if (is.null(channel_names <- names(observations))) {
        channel_names <- paste("Channel", seq_len(n_channels))
      }
    } else if (length(channel_names) != n_channels) {
      warning_("The length of {.arg channel_names} does not match the number 
               of channels. Names were not used.")
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

.check_transition_probs <- function(transition_probs, state_names = NULL) {
  stopifnot_(
    is.matrix(transition_probs),
    "{.arg transition_probs} is not a {.cls matrix}."
  )
  stopifnot_(
    dim(transition_probs)[1] == dim(transition_probs)[2],
    "{.arg transition_probs} is not a square {.cls matrix}."
  )
  n_states <- nrow(transition_probs)
  if (is.null(state_names)) {
    if (is.null(state_names <- rownames(transition_probs))) {
      state_names <- paste("State", seq_len(n_states))
    }
  } else {
    stopifnot_(
      length(state_names) == n_states,
      "Length of {.arg state_names} is not equal to the number of hidden 
      states."
    )
  }
  stopifnot_(
    isTRUE(all.equal(
      rowSums(transition_probs), 
      rep(1, dim(transition_probs)[1]), check.attributes = FALSE)),
    "Transition probabilities in {.arg transition_probs} do not sum to one."
  )
  dimnames(transition_probs) <- list(from = state_names, to = state_names)
  transition_probs
}

.check_initial_probs <- function(initial_probs, n_states, state_names = NULL) {
  stopifnot_(
    is.vector(initial_probs),
    "{.arg initial_probs} is not a {.cls vector}."
  )
  stopifnot_(
    length(initial_probs) == n_states,
    "Length of {.arg initial_probs} is not equal to the number of hidden 
      states."
  )
  
  stopifnot_(
    isTRUE(all.equal(sum(initial_probs), 1, check.attributes = FALSE)),
    "Initial state probabilities in {.arg initial_probs} do not sum to one."
  )
  names(initial_probs) <- state_names
  initial_probs
}

.check_emission_probs <- function(
    emission_probs, n_states, n_channels, n_symbols, state_names, symbol_names,
    channel_names = NULL) {
  
  if (n_channels == 1) {
    emission_probs <- list(emission_probs)
  }
  stopifnot_(
    length(emission_probs) == n_channels,
    "Number of channels defined by {.arg emission_probs} differs from one 
    defined by observations."
  )
  for (j in seq_len(n_channels)) {
    z <- if (n_channels > 1) paste0(" for channel ", j) else ""
    stopifnot_(
      is.matrix(emission_probs[[j]]),
      "{.arg emission_probs}{z} is not a {.cls matrix}." 
    )
    stopifnot_(
      isTRUE(all.equal(
        rowSums(emission_probs[[j]]), 
        rep(1, n_states), 
        check.attributes = FALSE)
      ),
      "{.arg emission_probs}{z} do not sum to one."
    )
    stopifnot_(
      nrow(emission_probs[[j]]) == n_states,
      "Number of rows in {.arg emission_probs}{z} is not equal to the number of 
      states."
    )
    stopifnot_(
      ncol(emission_probs[[j]]) == n_symbols[j],
      "Number of columns in {.arg emission_probs}{z} is not equal to the 
      number of symbols."
    )
  }
  
  if (is.null(channel_names)) {
    if (is.null(channel_names <- names(emission_probs))) {
      channel_names <- paste("Channel", 1:n_channels)
    }
  } else {
    if (length(channel_names) != n_channels) {
      warning_(
        "The length of {.arg channel_names} does not match the number of 
      channels. Names were not used."
      )
      channel_names <- paste("Channel", 1:n_channels)
    }
  }
  for (i in seq_len(n_channels)) {
    dimnames(emission_probs[[i]]) <- 
      list(state_names = state_names, symbol_names = symbol_names[[i]])
  }
  names(emission_probs) <- channel_names
  if (is.list(emission_probs) && length(emission_probs) == 1L) {
    emission_probs <- emission_probs[[1]]
  }
  emission_probs
}

.check_data <- function(data, id, time) {
  stopifnot_(
    is.data.frame(data), 
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(!missing(time), "Argument {.arg time} is missing.")
  stopifnot_(!missing(id), "Argument {.arg id} is missing.")
  stopifnot_(
    checkmate::test_string(x = time), 
    "Argument {.arg time} must be a single character string"
  )
  stopifnot_(
    checkmate::test_string(x = id), 
    "Argument {.arg id} must be a single character string"
  )
  stopifnot_(
    !is.null(data[[id]]), 
    "Can't find grouping variable {.var {id}} in {.arg data}."
  )
  stopifnot_(
    !is.null(data[[time]]), 
    "Can't find time index variable {.var {time}} in {.arg data}."
  )
  data <- data[order(data[[id]], data[[time]]), ]
  fill_time(data, id, time)
}