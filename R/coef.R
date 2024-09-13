#' Get the Estimated Regression Coefficients of Non-Homogeneous Hidden Markov 
#' Models
#' 
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`.
#' @param ... Ignored.
#' @rdname coef
#' @export
coef.nhmm <- function(object, probs = c(0.025, 0.975), ...) {
  S <- object$n_states
  M <- object$n_symbols
  beta_i_raw <- c(object$coefficients$beta_i_raw)
  beta_s_raw <- c(object$coefficients$beta_s_raw)
  beta_o_raw <- c(object$coefficients$beta_o_raw)
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
      observation = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = beta_o_raw
    )
  } else {
    beta_o <- data.frame(
      state = object$state_names,
      observation = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * (M[i] - 1))
      })),
      estimate = beta_o_raw
    )
  }
  
  stopifnot_(
    checkmate::test_numeric(
      x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
  )
  p_i <- length(beta_i_raw)
  p_s <- length(beta_s_raw)
  p_o <- length(beta_o_raw)
  sds <- try(
    diag(solve(-object$estimation_results$hessian)), 
    silent = TRUE
  )
  if (inherits(sds, "try-error")) {
    warning_(
      paste0(
        "Standard errors could not be computed due to singular Hessian. ",
        "Confidence intervals will not be provided."
      )
    )
    sds <- rep(NA, p_i + p_s + p_o)
  } else {
    if (any(sds < 0)) {
      warning_(
        paste0(
          "Standard errors could not be computed due to negative variances. ",
          "Confidence intervals will not be provided."
        )
      )
      sds <- rep(NA, p_i + p_s + p_o)
    } else {
      sds <- sqrt(sds)
    }
  }
  for(i in seq_along(probs)) {
    q <- qnorm(probs[i])
    beta_i[paste0("q", 100 * probs[i])] <- beta_i_raw + q * sds[seq_len(p_i)]
    beta_s[paste0("q", 100 * probs[i])] <- beta_s_raw + q * sds[p_i + seq_len(p_s)]
    beta_o[paste0("q", 100 * probs[i])] <- beta_o_raw + q * sds[p_i + p_s + seq_len(p_o)]
  }
  
  list(
    beta_initial = beta_i, 
    beta_transition = beta_s, 
    beta_emission = beta_o
  )
}
#' @rdname coef
#' @export
coef.mnhmm <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  
  S <- object$n_states
  M <- object$n_symbols
  D <- object$n_clusters
  beta_i_raw <- c(object$coefficients$beta_i_raw)
  beta_s_raw <- c(object$coefficients$beta_s_raw)
  beta_o_raw <- c(object$coefficients$beta_o_raw)
  theta_raw <- c(object$coefficients$theta_raw)
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
      observations = rep(object$symbol_names[-1], each = S * D),
      parameter = rep(object$coef_names_emission, each = S * (M - 1) * D),
      estimate = beta_o_raw,
      cluster = object$cluster_names
    )
  } else {
    beta_o <- data.frame(
      state = rep(object$state_names, each = D),
      observations = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S * D),
      parameter = rep(unlist(lapply(seq_len(object$n_channels), function(i) {
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
  
  stopifnot_(
    checkmate::test_numeric(
      x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
  )
  p_i <- length(beta_i_raw)
  p_s <- length(beta_s_raw)
  p_o <- length(beta_o_raw)
  p_c <- length(theta_raw)
  sds <- try(
    sqrt(diag(solve(-object$estimation_results$hessian))), 
    silent = TRUE
  )
  if (inherits(sds, "try-error")) {
    warning_(
      "Standard errors could not be computed due to singular Hessian. 
      Confidence intervals will not be provided."
    )
    sds <- rep(NA, p_i + p_s + p_o + p_c)
  }
  
  for(i in seq_along(probs)) {
    q <- qnorm(probs[i])
    beta_i[paste0("q", 100 * probs[i])] <- beta_i_raw + q * sds[seq_len(p_i)]
    beta_s[paste0("q", 100 * probs[i])] <- beta_s_raw + q * sds[p_i + seq_len(p_s)]
    beta_o[paste0("q", 100 * probs[i])] <- beta_o_raw + q * sds[p_i + p_s + seq_len(p_o)]
    theta[paste0("q", 100 * probs[i])] <- theta_raw + q * sds[p_i + p_s + p_o + seq_len(p_c)]
  }
  
  list(
    beta_initial = beta_i, 
    beta_transition = beta_s, 
    beta_emission = beta_o,
    beta_cluster = theta
  )
}
