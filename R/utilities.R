#' Quantile function with more data.frame friendly names
#' @noRd
quantileq <- function(x, probs, ...) {
  stats::setNames(
    stats::quantile(x, probs = probs, ...), 
    paste0("q", 100 * probs)
  )
}
#' Convert return code from estimate_nhmm and estimate_mnhmm to text
#' 
#' @param code Integer return code from `model$estimation_results$return_code`.
#' @return Code translated to informative message.
return_msg <- function(code) {
  
  msg <- paste0("Code ", code, "is not possible.")
  if (code < -1000) {
    x <- "Initial EM step failed\n. "
    code <- code + 1000
  } else {
    x <- ""
  }
  gamma <- NULL
  if (code %in% -c(100, 111:114)) gamma <- "gamma_pi"
  if (code %in% -c(200, 211:214)) gamma <- "gamma_A"
  if (code %in% -c(300, 311:314)) gamma <- "gamma_B"
  if (code %in% -c(400, 411:414)) gamma <- "gamma_omega"
  if (!is.null(gamma) && code %in% (-100 * (1:4))) {
    return(paste0(x,
                  "Error in M-step of ", gamma, " encountered expected count of zero. ",
                  "Try increasing the regularization via lambda or adjust parameter bounds ",
                  "to avoid extreme probabilities.")
    )
  }
  
  if (!(code %in% -(1:4))) {
    mstep <- paste0("Error in M-step of ", gamma, ". ")
  } else {
    mstep <- ""
  }
  
  e <- c(0, seq(110, 410, by = 100))
  if (code %in% -(1 + e)) {
    msg <- paste0(mstep, "NLOPT_FAILURE: Generic failure code.")
  }
  if (code %in% -(2 + e)) {
    msg <- paste0(
      mstep, 
      "NLOPT_INVALID_ARGS: Invalid arguments (e.g., lower bounds are ",
      "bigger than upper bounds, an unknown algorithm was specified)."
    )
  }
  if (code %in% -(3 + e)) {
    msg <- paste0(
      mstep, "NLOPT_OUT_OF_MEMORY: Ran out of memory."
    )
  }
  if (code %in% -(4 + e)) {
    msg <- paste0(
      mstep,
      "NLOPT_ROUNDOFF_LIMITED: Halted because roundoff errors limited progress."
    )
  }
  if (code == 0) {
    msg <- "Generic success."
  }
  if (code == 1) {
    msg <- "Generic success."
  }
  if (code == 2) {
    msg <- "Optimization stopped because stopval was reached."
  }
  if (code == 3) {
    msg <- "Optimization stopped because ftol_rel or ftol_abs was reached."
  }
  if (code == 4) {
    msg <- "Optimization stopped because xtol_rel or xtol_abs was reached."
  }
  if (code == 5) {
    msg <- "Optimization stopped because maxeval was reached."
  }
  if (code == 6) {
    msg <- " Optimization stopped because maxtime was reached."  
  }
  if (code == 7) {
    msg <- paste0(
      "NLopt terminated with generic error code -1, ", 
      "but the the relative or absolute change of the log-likelihood was ",
      "smaller than the tolerances so likely converged successfully."
    )
  }
  paste0(x, msg)
}

#' Split mnhmm for multiple nhmms for get_transitions etc
#' 
#' @noRd
split_mnhmm <- function(x) {
  D <- x$n_clusters
  models <- lapply(seq_len(D), \(i) {
    z <- x
    z$etas[1:3] <- lapply(x$etas[1:3], "[[", i)
    z$gammas[1:3] <- lapply(x$gammas[1:3], "[[", i)
    z$boot$gamma_pi <- lapply(x$boot$gamma_pi, "[[", i)
    z$boot$gamma_A <- lapply(x$boot$gamma_A, "[[", i)
    z$boot$gamma_B <- lapply(x$boot$gamma_B, "[[", i)
    z$state_names <- x$state_names[[i]]
    z$n_clusters <- 1L
    class(z) <- c(if(inherits(x, "fanhmm")) "fanhmm", "nhmm")
    z
  })
  names(models) <- x$cluster_names
  models
}

io <- function(X) {
  out <- attr(X, "icpt_only")
  if (is.null(out)) { 
    out <- vapply(X, io, TRUE)
  }
  out
}
iv <- function(X) {
  out <- attr(X, "iv")
  if (is.null(out)) { 
    out <- vapply(X, iv, TRUE)
  }
  out
}
tv <- function(X) {
  out <- attr(X, "tv")
  if (is.null(out)) { 
    out <- vapply(X, tv, TRUE)
  }
  out
}
K <- function(X) {
  out <- length(attr(X, "coef_names"))
  if (out == 0) { 
    out <- vapply(X, K, 1L)
  }
  out
}
#' Check if x is a list or a list of length n consisting of lists
#' @noRd
is_list_of_lists <- function(x, n) {
  is.list(x) && length(x) == n && all(vapply(x, is.list, TRUE))
}
#' Check that x is stslist or list of length n consisting of stslist objects
#' @noRd
is_stslist <- function(x, n) {
  if (TraMineR::is.stslist(x)) {
    return(n == 1L)
  } else if (is.list(x)) {
    return(all(vapply(x, \(i) TraMineR::is.stslist(i), TRUE)) && length(x) == n)
  }
  FALSE
}
#' Check that x is formula or list of length n consisting of formulas 
#' @noRd
is_formula <- function(x, n) {
  if (inherits(x, "formula")) {
    return(n == 1L)
  } else if (is.list(x)) {
    return(all(vapply(x, \(i) inherits(i, "formula"), TRUE)) && length(x) == n)
  }
  FALSE
}
#' check that x is valid probability mass function
#' @noRd
is_pmf <- function(x, n, tol = 1e-8) {
  is.numeric(x) && length(x) == n && all(is.finite(x)) && all(x >= 0) && abs(sum(x) - 1) < tol
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
create_obs <- function(model) {
  # avoid CRAN check warnings due to NSE
  cols <- NULL
  id <- model$data[[model$id_variable]]
  obs <- qM(model$data[, cols, env = list(cols = I(model$responses))]) - 1L
  for (i in seq_len(ncol(obs))) {
    obs[is.na(obs[, i]), i] <- model$n_symbols[i]
  }
  lapply(rsplit(obs, fl = id), \(y) t(y))
}
#' Create obsArray for various pre-2.0.0 C++ functions
#' 
#' @noRd
create_obsArray <- function(model) {
  obsArray <- array(
    0L, 
    c(model$n_sequences, model$length_of_sequences, model$n_channels)
  )
  if (model$n_channels == 1) {
    model$observations <- list(model$observations)
  }
  for (i in seq_len(model$n_channels)) {
    obsArray[, , i] <- unlist(lapply(model$observations[[i]], as.integer)) - 1L
    obsArray[, , i][obsArray[, , i] > model$n_symbols[i]] <- model$n_symbols[i]
  }
  aperm(obsArray)
}
#' Create emissionArray for various pre-2.0.0 C++ functions
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

#' Count unique values
#' @noRd
n_unique <- function(x) {length(unique(x))}

#' Create a factor from character vector of levels in order of appearance
#' @noRd
as_factor <- function(x) {
  if (is.factor(x)) {
    x
  } else {
    factor(x, levels = unique(x))
  }
}
