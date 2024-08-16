#' Stop Function Execution Unless Condition Is True
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See [cli::cli_abort()].
#' @noRd
stopifnot_ <- function (cond, message, ..., call = rlang::caller_env()) {
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
stop_ <- function (message, ..., call = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
}
#' Checks Whether the Formula Contains Only the Intercept Term
#' 
#' @param f A formula object.
#' @noRd
intercept_only <- function(f) {
  identical(deparse(update(f, 0 ~ .)), "0 ~ 1")
}
#' @noRd
check_positive_integer <- function(x, name = "threads") {
  stopifnot_(
    x > 0 && x == as.integer(x),
    "{.arg {name}} must be a positive integer."
  )
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