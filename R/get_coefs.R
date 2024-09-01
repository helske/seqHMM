#' Get the Estimated Regression Coefficients of Non-Homogeneous Hidden Markov 
#' Models
#' 
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param nsim Non-negative integer defining the number of samples from the 
#' normal approximation of the model parameters used in 
#' computing the approximate quantiles of the estimates. If `0`, only point 
#' estimates are returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.975)`.
#' @rdname coef
#' @export
coef.nhmm <- function(object, nsim = 0, probs = c(0.025, 0.975), ...) {
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  S <- object$n_states
  M <- object$n_symbols
  beta_i_raw <- c(object$estimation_results$parameters$beta_i_raw)
  beta_s_raw <- c(object$estimation_results$parameters$beta_s_raw)
  beta_o_raw <- c(object$estimation_results$parameters$beta_o_raw)
  beta_i <- data.frame(
    state = object$state_names[-1],
    parameter = rep(object$coef_names_initial, each = S - 1),
    estimate = beta_i_raw
  )
  beta_s <- data.frame(
    state_from = object$state_names,
    state_to = rep(object$state_names[-1], each = S),
    parameter = rep(object$coef_names_transition, each = S * (S - 1)),
    estimate = beta_s_raw
  )
  if (object$n_channels == 1) {
    beta_o <- data.frame(
      state = object$state_names,
      symbol = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = beta_o_raw
    )
  } else {
    beta_o <- data.frame(
      state = object$state_names,
      symbol = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S),
      parameter = unlist(lapply(seq_len(model$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * (M[i] - 1))
      })),
      estimate = beta_o_raw
    )
  }
  if (nsim > 0) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    chol_precision <- chol(-object$estimation$hessian)
    U <- backsolve(chol_precision, diag(ncol(chol_precision)))
    x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
    x <- t(sweep(x, 2, c(beta_i_raw, beta_s_raw, beta_o_raw), "+"))
    p_i <- length(beta_i_raw)
    p_s <- length(beta_s_raw)
    p_o <- length(beta_o_raw)
    quantiles <- fast_quantiles(x[seq_len(p_i), ], probs)
    for(i in seq_along(probs)) {
      beta_i[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    quantiles <- fast_quantiles(x[p_i + seq_len(p_s), ], probs)
    for(i in seq_along(probs)) {
      beta_s[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    quantiles <- fast_quantiles(x[p_i + p_s + seq_len(p_o), ], probs)
    for(i in seq_along(probs)) {
      beta_o[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
  }
  list(
    beta_initial = beta_i, 
    beta_transition = beta_s, 
    beta_emission = beta_o
  )
}
#' @rdname coef
#' @export
coef.mnhmm <- function(object, nsim = 0, probs = c(0.025, 0.975), ...) {
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  S <- object$n_states
  M <- object$n_symbols
  D <- object$n_clusters
  beta_i_raw <- c(object$estimation_results$parameters$beta_i_raw)
  beta_s_raw <- c(object$estimation_results$parameters$beta_s_raw)
  beta_o_raw <- c(object$estimation_results$parameters$beta_o_raw)
  theta_raw <- c(object$estimation_results$parameters$theta_raw)
  beta_i <- data.frame(
    state = rep(object$state_names[-1], each = D),
    parameter = rep(object$coef_names_initial, each = (S - 1) * D),
    estimate = beta_i_raw,
    cluster = object$cluster_names
  )
  beta_s <- data.frame(
    state_from = rep(object$state_names, each = D),
    state_to = rep(object$state_names[-1], each = S * D),
    parameter = rep(object$coef_names_transition, each = S * (S - 1) * D),
    estimate = beta_s_raw,
    cluster = object$cluster_names
  )
  if (object$n_channels == 1) {
    beta_o <- data.frame(
      state = rep(object$state_names, each = D),
      symbol = rep(object$symbol_names[-1], each = S * D),
      parameter = rep(object$coef_names_emission, each = S * (M - 1) * D),
      estimate = beta_o_raw,
      cluster = object$cluster_names
    )
  } else {
    beta_o <- data.frame(
      state = rep(object$state_names, each = D),
      symbol = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S * D),
      parameter = rep(unlist(lapply(seq_len(model$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * (M[i] - 1))
      })), each = D),
      estimate = beta_o_raw
    )
  }
  theta <- data.frame(
    cluster = object$cluster_names[-1],
    parameter = rep(object$coef_names_cluster, each = D - 1),
    estimate = theta_raw
  )
  if (nsim > 0) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    chol_precision <- chol(-object$estimation$hessian)
    U <- backsolve(chol_precision, diag(ncol(chol_precision)))
    x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
    x <- t(sweep(x, 2, c(beta_i_raw, beta_s_raw, beta_o_raw, theta_raw), "+"))
    p_i <- length(beta_i_raw)
    p_s <- length(beta_s_raw)
    p_o <- length(beta_o_raw)
    p_c <- length(theta_o_raw)
    quantiles <- fast_quantiles(x[seq_len(p_i), ], probs)
    for(i in seq_along(probs)) {
      beta_i[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    quantiles <- fast_quantiles(x[p_i + seq_len(p_s), ], probs)
    for(i in seq_along(probs)) {
      beta_s[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    quantiles <- fast_quantiles(x[p_i + p_s + seq_len(p_o), ], probs)
    for(i in seq_along(probs)) {
      beta_o[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    quantiles <- fast_quantiles(x[p_i + p_s + p_o + seq_len(p_c), ], probs)
    for(i in seq_along(probs)) {
      theta[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
  }
  list(
    beta_initial = beta_i, 
    beta_transition = beta_s, 
    beta_emission = beta_o,
    beta_cluster = theta
  )
}