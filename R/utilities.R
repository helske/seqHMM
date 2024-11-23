no_itcp_idx <- function(x) {
  which(arrayInd(seq_along(x), dim(x))[ , 2] != 1)
}
#' Split mnhmm for multiple nhmms for get_probs
#' 
#' @noRd
split_mnhmm <- function(x) {
  D <- x$n_clusters
  models <- lapply(seq_len(D), function(i) {
    z <- x
    z$etas[1:3] <- lapply(x$etas[1:3], "[[", i)
    z$gammas[1:3] <- lapply(x$gammas[1:3], "[[", i)
    z$boot$gamma_pi <- lapply(x$boot$gamma_pi, "[[", i)
    z$boot$gamma_A <- lapply(x$boot$gamma_A, "[[", i)
    z$boot$gamma_B <- lapply(x$boot$gamma_B, "[[", i)
    z$state_names <- x$state_names[[i]]
    z$n_clusters <- 1L
    class(z) <- "nhmm"
    z
  })
  names(models) <- x$cluster_names
  models
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
#' (Regularized) Inverse of softmax(Q*eta)
#' 
#' @noRd
p_to_eta <-function(x) {
  Q <- create_Q(length(x))
  x <- pmin(pmax(x, 1e-6), 1-1e-6)
  x <- x / sum(x)
  log_x <- log(x)
  t(Q) %*% (log_x - mean(log_x))
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
#' Convert error message to text
#' @noRd
error_msg <- function(error) {
  
  if (error %in% c(1, -(101:104))) gamma <- "gamma_pi"
  if (error %in% c(2, -(201:204))) gamma <- "gamma_A"
  if (error %in% c(3, -(301:304))) gamma <- "gamma_B"
  if (error %in% c(4, -(401:404))) gamma <- "gamma_omega"
  
  nonfinite_msg <- paste0(
    "Error: Some of the values in ", gamma, " are nonfinite, likely due to ",
    "zero expected counts. Try increasing the penalty lambda or ",
    "pseudocounts (in case of EM algorithm) to avoid extreme probabilities.")
  if (!(error %in% -(1:4))) {
    mstep <- paste0("Error in M-step of ", gamma, ". ")
  } else {
    mstep <- ""
  }
  
  e <- seq(0, 400, by = 100)
  if (error %in% 1:4) {
    msg <- nonfinite_msg
  }
  if (error %in% (-1 - e)) {
    msg <- paste0(mstep, "NLOPT_FAILURE: Generic failure code.")
  }
  if (error %in% (-2 - e)) {
    msg <- paste0(
      mstep, 
      "NLOPT_INVALID_ARGS: Invalid arguments (e.g., lower bounds are ",
      "bigger than upper bounds, an unknown algorithm was specified)."
    )
  }
  if (error %in% (-3 - e)) {
    msg <- paste0(
      mstep, "NLOPT_OUT_OF_MEMORY: Ran out of memory."
    )
  }
  if (error %in% (-4 - e)) {
    msg <- paste0(
      mstep,
      "NLOPT_ROUNDOFF_LIMITED: Halted because roundoff errors limited progress."
    )
  }
  msg
}