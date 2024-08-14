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

#' Checks Whether the Formula Contains Only the Intercept Term
#' 
#' @param f A formula object.
#' @noRd
intercept_only <- function(f) {
  identical(deparse(update(f, 0 ~ .)), "0 ~ 1")
}
