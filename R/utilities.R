#' Filter data.frame Rows with Time not Exceeding Sequence Length
#' @noRd
non_void_rows <- function(d, time, id, sequence_lengths) {
  d[d[[time]] <= sequence_lengths[match(d[[id]], unique(d[[id]]))], ]
}
#' Check If an Object Is a List of Lists
#' 
#' @noRd
is_list_of_lists <- function(x) {
  if (!is.list(x)) {
    return(FALSE)
  } else {
    if (all(unlist(lapply(x, is.list)))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
#' Regularized Inverse Softmax Function
#' 
#' @noRd
inv_softmax <- function(x) {
  x <- pmin(pmax(x, 0.001), 0.999)
  x <- x / sum(x)
  log(x) - log(x[1])
}

#' Stop Function Execution Unless Condition Is True
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See [cli::cli_abort()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#' 
stopifnot_ <- function(cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}
#' Give a Informational Message
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_inform()].
#' @param ... See [cli::cli_inform()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
message_ <- function(message, ...) {
  cli::cli_inform(message, ..., .envir = parent.frame())
}
#' Issue a Warning
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_warn()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#' 
warning_ <- function (message, ...) {
  cli::cli_warn(message, ..., .envir = parent.frame())
}
#' Issue an Error
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See [cli::cli_abort()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#' 
stop_ <- function(message, ..., call = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
}
#' Checks Whether the Formula Contains Only the Intercept Term
#' 
#' @param f A formula object.
#' @noRd
intercept_only <- function(f) {
  identical(deparse(update(f, 0 ~ .)), "0 ~ 1")
}
#' Create obsArray for Various C++ functions
#' 
#' @noRd
create_obsArray <- function(model) {
  obsArray <- array(
    0L, 
    c(model$n_sequences, model$length_of_sequences, model$n_channels))
  if (model$n_channels == 1) {
    model$observations <- list(model$observations)
  }
  for (i in 1:model$n_channels) {
    obsArray[, , i] <- unlist(lapply(model$observations[[i]], as.integer)) - 1L
    obsArray[, , i][obsArray[, , i] > model$n_symbols[i]] <- model$n_symbols[i]
    stopifnot_(
      sum(obsArray[, , i] < model$n_symbols[i]) > 0,
      "One channel contains only missing values, model is degenerate."
    )
  }
  aperm(obsArray)
}
#' Create emissionArray for Various C++ functions
#' 
#' @noRd
create_emissionArray <- function(model) {
  emissionArray <- array(1, c(model$n_states, max(model$n_symbols) + 1, model$n_channels))
  if (model$n_channels == 1) {
    model$emission_probs <- list(model$emission_probs)
  }
  for (i in 1:model$n_channels) {
    emissionArray[, 1:model$n_symbols[i], i] <- model$emission_probs[[i]]
  }
  emissionArray
}

#' Remove rows corresponding to void symbols
#' @noRd
remove_voids <- function(model, x) {
  
  x_id <- x[[model$id_variable]]
  x_time <- x[[model$time_variable]]
  if (model$n_channels == 1) {
    time <- colnames(model$observations)
    id <- rownames(model$observations)
  } else {
    time <- colnames(model$observations[[1]])
    id <- rownames(model$observations[[1]])
  }
  do.call(
    rbind, 
    lapply(seq_len(model$n_sequences), function(i) {
      idx <- which(
        x_id == id[i] & x_time %in% time[seq_len(model$sequence_lengths[i])]
      )
      x[idx, ]
    })
  )
}
