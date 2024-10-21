#' Create observations for the model objects
#' 
#' Note that for backward compatibility reasons `length_of_sequences` refers 
#' to the maximum length of sequences, whereas `sequence_lengths` refers to 
#' the actual non-void lengths of each sequence.
#'@noRd
#'
.check_observations <- function(x, channel_names = NULL, 
                                nhmm = FALSE) {
  
  if (TraMineR::is.stslist(x)) {
    multichannel <- FALSE
  } else {
    if (nhmm) {
      stopifnot_(
        is_list_of_lists(x) && all(unlist(lapply(x, TraMineR::is.stslist))),
        "{.arg observations} should be a {.cls stslist} object created with 
        {.fn seqdef}, a {.cls list} of {.cls stslist} objects, or a 
        {.cls character} vector containing names of the response variables in 
        {.arg data}."
      )
    } else {
      stopifnot_(
        is_list_of_lists(x) && all(unlist(lapply(x, TraMineR::is.stslist))),
        "{.arg observations} should be a {.cls stslist} object created with 
        {.fn seqdef}, or a {.cls list} of {.cls stslist} objects in a 
        multichannel case."
      )
    }
    multichannel <- TRUE
    stopifnot_(
      length(unique(sapply(x, nrow))) == 1,
      "The number of subjects (rows) is not the same in all channels."
    )
    stopifnot_(
      length(unique(sapply(x, ncol))) == 1,
      "The length of the sequences (columns) is not the same in all channels."
    )
  }
  # Single channel but observations is a list
  if (is.list(x) && 
      !TraMineR::is.stslist(x) && 
      length(x) == 1) {
    x <- x[[1]]
    multichannel <- FALSE
  }
  n_channels <- ifelse(multichannel, length(x), 1L)
  if (multichannel) {
    n_sequences <- nrow(x[[1]])
    length_of_sequences <- ncol(x[[1]])
    
    symbol_names <- lapply(x, alphabet)
    n_symbols <- lengths(symbol_names)
    dim(n_symbols) <- length(n_symbols)
    if (is.null(channel_names)) {
      if (is.null(channel_names <- names(x))) {
        channel_names <- paste("Channel", seq_len(n_channels))
      }
    } else if (length(channel_names) != n_channels) {
      warning_("The length of {.arg channel_names} does not match the number 
               of channels. Names were not used.")
      channel_names <- paste("Channel", seq_len(n_channels))
    }
    sequence_lengths <- do.call(pmax, lapply(x, TraMineR::seqlength))
    times <- colnames(x[[1]])
    na_times <- suppressWarnings(any(is.na(timenames <- as.numeric(times))))
    if (na_times) {
      na_times <- suppressWarnings(any(is.na(timenames <- as.numeric(sub('.', '', times)))))
      if (na_times) {
        warning_(
          paste0(
            "Time indices (column names) of sequences are not coarceable to ",
            "numeric. Replacing them with integers."
          )
        )
        timenames <- seq_len(ncol(x[[1]]))
      }
    }
    stopifnot_(
      identical(sort(timenames), timenames),
      paste0(
        "The numeric time indices based on column names of sequence object ", 
        "are not numerically sorted. Please recode the column names.")
    )
    for (i in seq_len(length(x))) {
      colnames(x[[i]]) <- timenames 
    }
  } else {
    n_sequences <- nrow(x)
    length_of_sequences <- ncol(x)
    
    symbol_names <- alphabet(x)
    n_symbols <- length(symbol_names)
    
    if (is.null(channel_names)) {
      channel_names <- "Observations"
    }
    sequence_lengths <- TraMineR::seqlength(x)
    times <- colnames(x)
    na_times <- suppressWarnings(any(is.na(timenames <- as.numeric(times))))
    if (na_times) {
      na_times <- suppressWarnings(any(is.na(timenames <- as.numeric(sub('.', '', times)))))
      if (na_times) {
        warning_(
          paste0(
            "Time indices (column names) of sequences are not coarceable to ",
            "numeric. Replacing them with integers."
          )
        )
        timenames <- seq_len(ncol(x))
      }
    }
    stopifnot_(
      identical(sort(timenames), timenames),
      paste0(
        "The numeric time indices based on column names of sequence object ", 
        "are not numerically sorted. Please recode the column names.")
    )
    colnames(x) <- timenames 
  }
  sequence_lengths <- as.integer(sequence_lengths)
  dim(sequence_lengths) <- length(sequence_lengths)
  if (n_channels > 1L) {
    nobs <- sum(sapply(x, function(x) {
      sum(!(x == attr(x, "nr") | x == attr(x, "void") | is.na(x)))
    })) / n_channels
  } else {
    nobs <- sum(!(x == attr(x, "nr") | x == attr(x, "void") | is.na(x)))
  }
  
  attr(x, "n_channels") <- ifelse(
    multichannel, length(x), 1L)
  attr(x, "n_sequences") <- n_sequences
  attr(x, "length_of_sequences") <- length_of_sequences
  attr(x, "sequence_lengths") <- sequence_lengths
  attr(x, "symbol_names") <- symbol_names
  attr(x, "n_symbols") <- n_symbols
  attr(x, "channel_names") <- channel_names
  attr(x, "nobs") <- nobs
  x
}

.check_transition_probs <- function(x, state_names = NULL) {
  stopifnot_(
    is.matrix(x),
    "{.arg x} is not a {.cls matrix}."
  )
  stopifnot_(
    dim(x)[1] == dim(x)[2],
    "{.arg transition_probs} is not a square {.cls matrix}."
  )
  n_states <- nrow(x)
  if (is.null(state_names)) {
    if (is.null(state_names <- rownames(x))) {
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
      rowSums(x), 
      rep(1, dim(x)[1]), check.attributes = FALSE)),
    "Transition probabilities in {.arg transition_probs} do not sum to one."
  )
  dimnames(x) <- list(from = state_names, to = state_names)
  x
}

.check_initial_probs <- function(x, n_states, state_names = NULL) {
  stopifnot_(
    is.vector(x),
    "{.arg initial_probs} is not a {.cls vector}."
  )
  stopifnot_(
    length(x) == n_states,
    "Length of {.arg initial_probs} is not equal to the number of hidden 
      states."
  )
  
  stopifnot_(
    isTRUE(all.equal(sum(x), 1, check.attributes = FALSE)),
    "Initial state probabilities in {.arg x} do not sum to one."
  )
  names(x) <- state_names
  x
}

.check_emission_probs <- function(
    x, n_states, n_channels, n_symbols, state_names, symbol_names,
    channel_names = NULL) {
  
  if (n_channels == 1) {
    x <- list(x)
    symbol_names <- list(symbol_names)
  }
  stopifnot_(
    length(x) == n_channels,
    "Number of channels defined by {.arg emission_probs} differs from one 
    defined by observations."
  )
  for (j in seq_len(n_channels)) {
    z <- if (n_channels > 1) paste0(" for channel ", j) else ""
    stopifnot_(
      is.matrix(x[[j]]),
      "{.arg emission_probs}{z} is not a {.cls matrix}." 
    )
    stopifnot_(
      isTRUE(all.equal(
        rowSums(x[[j]]), 
        rep(1, n_states), 
        check.attributes = FALSE)
      ),
      "{.arg emission_probs}{z} do not sum to one."
    )
    stopifnot_(
      nrow(x[[j]]) == n_states,
      "Number of rows in {.arg emission_probs}{z} is not equal to the number of 
      states."
    )
    stopifnot_(
      ncol(x[[j]]) == n_symbols[j],
      "Number of columns in {.arg emission_probs}{z} is not equal to the 
      number of symbols."
    )
  }
  
  if (is.null(channel_names)) {
    if (is.null(channel_names <- names(x))) {
      channel_names <- paste("Channel", seq_len(n_channels))
    }
  } else {
    if (length(channel_names) != n_channels) {
      warning_(
        "The length of {.arg channel_names} does not match the number of 
      channels. Names were not used."
      )
      channel_names <- paste("Channel", seq_len(n_channels))
    }
  }
  for (i in seq_len(n_channels)) {
    dimnames(x[[i]]) <- 
      list(state_names = state_names, symbol_names = symbol_names[[i]])
  }
  names(x) <- channel_names
  if (is.list(x) && length(x) == 1L) {
    x <- x[[1]]
  }
  x
}

.check_data <- function(data, time, id) {
  stopifnot_(
    is.data.frame(data), 
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(!missing(time), "Argument {.arg time} is missing.")
  stopifnot_(!missing(id), "Argument {.arg id} is missing.")
  stopifnot_(
    checkmate::test_string(x = time), 
    "Argument {.arg time} must be a single character string."
  )
  stopifnot_(
    checkmate::test_string(x = id), 
    "Argument {.arg id} must be a single character string."
  )
  stopifnot_(
    !is.null(data[[id]]), 
    "Can't find grouping variable {.var {id}} in {.arg data}."
  )
  stopifnot_(
    !is.null(data[[time]]), 
    "Can't find time index variable {.var {time}} in {.arg data}."
  )
  stopifnot_(
    is.numeric(data[[time]]),
    c(
      "Time index variable {.arg {time}} must be of type {.cls numeric} or 
      {.cls integer}."
    )
  )
  data <- data[order(data[[id]], data[[time]]), ]
  fill_time(data, id, time)
}
