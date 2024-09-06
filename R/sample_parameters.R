#' Samples parameters of non-homegeneous Markov hidden models using the normal approximation.
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param nsim A non-negative integer defining the number of samples from the 
#' normal approximation.
#' @param probs A numeric vector defining the quantiles of interest.
#' @param return_samples A logical indicating whether to return the samples 
#' instead of quantiles.
#' @noRd
sample_parameters <- function(model, nsim, probs, return_samples = FALSE) {
  stopifnot_(
    checkmate::test_int(x = nsim, lower = 1L), 
    "Argument {.arg nsim} must be a single positive integer."
  )
  stopifnot_(
    !return_samples && checkmate::test_numeric(
      x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
  )
  mixture <- inherits(model, "mnhmm")
  X_initial <- t(model$X_initial)
  X_transition <- aperm(model$X_transition, c(3, 1, 2))
  X_emission <- aperm(model$X_emission, c(3, 1, 2))
  if (mixture) X_cluster <- t(model$X_cluster)
  
  chol_precision <- chol(-model$estimation$hessian)
  U <- backsolve(chol_precision, diag(ncol(chol_precision)))
  x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
  beta_i_raw <- model$coefficientsbeta_i_raw
  beta_s_raw <- model$coefficientsbeta_s_raw
  beta_o_raw <- model$coefficientsbeta_o_raw
  pars <- c(beta_i_raw, beta_s_raw, beta_o_raw)
  p_i <- length(beta_i_raw)
  p_s <- length(beta_s_raw)
  p_o <- length(beta_o_raw)
  if (mixture) {
    theta_raw <- model$coefficientstheta_raw
    pars <- c(pars, theta_raw)
    if (mixture) p_d <- length(theta_raw)
  }
  x <- t(sweep(x, 2, pars, "+"))
  samples_pi <- apply(
    x[seq_len(p_i), ], 2, function(z) {
      z <- array(z, dim = dim(beta_i_raw))
      get_pi(z, X_initial)
    }
  )
  samples_A <- apply(
    x[p_i + seq_len(p_s), ], 2, function(z) {
      z <- array(z, dim = dim(beta_s_raw))
      unlist(get_A(aperm(z, c(2, 3, 1)), X_transition))
    }
  )
  samples_B <- apply(
    x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
      z <- array(z, dim = dim(beta_o_raw))
      unlist(get_B(aperm(z, c(2, 3, 1)), X_emission))
    }
  )
  if (mixture) {
    samples_omega <- apply(
      x[p_i + p_s + p_o + seq_len(p_d), ], 2, function(z) {
        z <- array(z, dim = dim(theta_raw))
        unlist(get_omega(z, X_cluster))
      }
    )
  }
  if (return_samples) {
    list(
      samples_pi = samples_pi,
      samples_A = samples_A,
      samples_B = samples_B,
      samples_omega = if (mixture) samples_omega
    )
  } else {
    list(
      quantiles_pi = fast_quantiles(samples_pi, probs),
      quantiles_A = fast_quantiles(samples_A, probs),
      quantiles_B = fast_quantiles(samples_B, probs),
      quantiles_omega = if (mixture) fast_quantiles(samples_omega, probs)
    )
  }
}
