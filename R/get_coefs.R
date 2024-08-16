#' Get the Estimated Regression Coefficients NHMM
#' 
#' @param object An object of class `nhmm`.
#' @param nsim Non-negative integer defining the number of samples used in 
#' constructing confidence intervals. If `0`, only point estimates are 
#' returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.975)`.
#' @export
coef.nhmm <- function(object, nsim = 0, probs = c(0.025, 0.975)) {
  
  S <- object$n_states
  M <- object$n_symbols
  beta_i <- data.frame(
    type = "Initial probability",
    state = object$state_names[-1],
    parameter = rep(object$coef_names_initial, each = S - 1),
    estimate = c(object$estimation_results$parameters$beta_i_raw)
  )
  beta_s <- data.frame(
    type = "Transition probability",
    state_from = object$state_names,
    state_to = rep(object$state_names[-1], each = S),
    parameter = rep(object$coef_names_transition, each = S * (S - 1)),
    estimate = c(object$estimation_results$parameters$beta_s_raw)
  )
  if (model$n_channels == 1) {
    beta_o <- data.frame(
      type = "Emission probability",
      state = object$state_names,
      symbol = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = c(object$estimation_results$parameters$beta_o_raw)
    )
  } else {
    beta_o <- data.frame(
      type = "Emission probability",
      state = object$state_names,
      symbol = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = c(object$estimation_results$parameters$beta_o_raw)
    )
  }
  if (nsim > 0) {
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
  rbind(beta_i, beta_s, beta_o)
}