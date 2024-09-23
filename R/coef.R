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
  gamma_pi_raw <- c(object$coefficients$gamma_pi_raw)
  gamma_pi <- data.frame(
    state = object$state_names[-1],
    parameter = rep(object$coef_names_initial, each = S - 1),
    estimate = gamma_pi_raw
  )
  gamma_A_raw <- c(object$coefficients$gamma_A_raw)
  gamma_A <- data.frame(
    state_from = object$state_names,
    state_to = rep(object$state_names[-1], each = S),
    parameter = rep(object$coef_names_transition, each = S * (S - 1)),
    estimate = gamma_A_raw
  )
  if (object$n_channels == 1) {
    gamma_B_raw <- c(object$coefficients$gamma_B_raw)
    gamma_B <- data.frame(
      state = object$state_names,
      observation = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = gamma_B_raw
    )
  } else {
    gamma_B_raw <- unlist(object$coefficients$gamma_B_raw)
    gamma_B <- data.frame(
      state = object$state_names,
      observation = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * (M[i] - 1))
      })),
      estimate = gamma_B_raw
    )
  }
  
  # stopifnot_(
  #   checkmate::test_numeric(
  #     x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
  #   ),
  #   "Argument {.arg probs} must be a {.cls numeric} vector with values
  #     between 0 and 1."
  # )
  # p_i <- length(gamma_pi_raw)
  # p_s <- length(gamma_A_raw)
  # p_o <- length(gamma_B_raw)
  # sds <- try(
  #   diag(solve(-object$estimation_results$hessian)), 
  #   silent = TRUE
  # )
  # if (inherits(sds, "try-error")) {
  #   warning_(
  #     paste0(
  #       "Standard errors could not be computed due to singular Hessian. ",
  #       "Confidence intervals will not be provided."
  #     )
  #   )
  #   sds <- rep(NA, p_i + p_s + p_o)
  # } else {
  #   if (any(sds < 0)) {
  #     warning_(
  #       paste0(
  #         "Standard errors could not be computed due to negative variances. ",
  #         "Confidence intervals will not be provided."
  #       )
  #     )
  #     sds <- rep(NA, p_i + p_s + p_o)
  #   } else {
  #     sds <- sqrt(sds)
  #   }
  # }
  # for(i in seq_along(probs)) {
  #   q <- qnorm(probs[i])
  #   gamma_pi[paste0("q", 100 * probs[i])] <- gamma_pi_raw + q * sds[seq_len(p_i)]
  #   gamma_A[paste0("q", 100 * probs[i])] <- gamma_A_raw + q * sds[p_i + seq_len(p_s)]
  #   gamma_B[paste0("q", 100 * probs[i])] <- gamma_B_raw + q * sds[p_i + p_s + seq_len(p_o)]
  # }
  
  list(
    gamma_pinitial = gamma_pi, 
    beta_transition = gamma_A, 
    beta_emission = gamma_B
  )
}
#' @rdname coef
#' @export
coef.mnhmm <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  
  S <- object$n_states
  M <- object$n_symbols
  D <- object$n_clusters
  gamma_pi_raw <- unlist(object$coefficients$gamma_pi_raw)
  K_i <- length(object$coef_names_initial)
  gamma_pi <- data.frame(
    state = unlist(lapply(object$state_names, function(x) x[-1])),
    parameter = rep(object$coef_names_initial, each = (S - 1)),
    estimate = unlist(gamma_pi_raw),
    cluster = rep(object$cluster_names, each = (S - 1) * K_i)
  )
  gamma_A_raw <- unlist(object$coefficients$gamma_A_raw)
  K_s <- length(object$coef_names_transition)
  gamma_A <- data.frame(
    state_from = unlist(object$state_names),
    state_to = rep(
      unlist(lapply(object$state_names, function(x) x[-1])), 
      each = S
    ),
    parameter = rep(object$coef_names_transition, each = S * (S - 1)),
    estimate = unlist(gamma_A_raw),
    cluster = rep(object$cluster_names, each = (S - 1) * S * K_s)
  )
  K_o <- length(object$coef_names_emission)
  if (object$n_channels == 1) {
    gamma_B <- data.frame(
      state = unlist(object$state_names),
      observations = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = unlist(gamma_B_raw),
      cluster =  rep(object$cluster_names, each = S * (S - 1) * K_o)
    )
  } else {
    gamma_B <- data.frame(
      state = unlist(object$state_names),
      observations = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * (M[i] - 1))
      })),
      estimate = unlist(gamma_B_raw),
      cluster =  unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$cluster_names, each = S * (M[i] - 1) * K_o)
      }))
    )
  }
  gamma_omega_raw <- c(object$coefficients$gamma_omega_raw)
  gamma_omega <- data.frame(
    cluster = object$cluster_names[-1],
    parameter = rep(object$coef_names_cluster, each = D - 1),
    estimate = gamma_omega_raw
  )
  
  # stopifnot_(
  #   checkmate::test_numeric(
  #     x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
  #   ),
  #   "Argument {.arg probs} must be a {.cls numeric} vector with values
  #     between 0 and 1."
  # )
  # p_i <- length(gamma_pi_raw)
  # p_s <- length(gamma_A_raw)
  # p_o <- length(gamma_B_raw)
  # p_d <- length(gamma_omega_raw)
  # sds <- try(
  #   sqrt(diag(solve(-object$estimation_results$hessian))), 
  #   silent = TRUE
  # )
  # if (inherits(sds, "try-error")) {
  #   warning_(
  #     "Standard errors could not be computed due to singular Hessian. 
  #     Confidence intervals will not be provided."
  #   )
  #   sds <- rep(NA, p_i + p_s + p_o + p_d)
  # }
  # 
  # for(i in seq_along(probs)) {
  #   q <- qnorm(probs[i])
  #   gamma_pi[paste0("q", 100 * probs[i])] <- gamma_pi_raw + q * sds[seq_len(p_i)]
  #   gamma_A[paste0("q", 100 * probs[i])] <- gamma_A_raw + q * sds[p_i + seq_len(p_s)]
  #   gamma_B[paste0("q", 100 * probs[i])] <- gamma_B_raw + q * sds[p_i + p_s + seq_len(p_o)]
  #   gamma_omega[paste0("q", 100 * probs[i])] <- gamma_omega_raw + q * sds[p_i + p_s + p_o + seq_len(p_d)]
  # }
  
  list(
    gamma_initial = gamma_pi, 
    gamma_transition = gamma_A, 
    gamma_emission = gamma_B,
    gamma_cluster = gamma_omega
  )
}
