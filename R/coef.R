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
  gamma_pi <- c(eta_to_gamma(object$coefficients$eta_pi))
  gamma_pi <- data.frame(
    state = object$state_names,
    parameter = rep(object$coef_names_initial, each = S),
    estimate = gamma_pi
  )
  eta_A <- c(object$coefficients$eta_A)
  gamma_A <- data.frame(
    state_from = object$state_names,
    state_to = rep(object$state_names[-1], each = S),
    parameter = rep(object$coef_names_transition, each = S * (S - 1)),
    estimate = eta_A
  )
  if (object$n_channels == 1) {
    eta_B <- c(object$coefficients$eta_B)
    gamma_B <- data.frame(
      state = object$state_names,
      observation = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = eta_B
    )
  } else {
    eta_B <- unlist(object$coefficients$eta_B)
    gamma_B <- data.frame(
      state = object$state_names,
      observation = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * (M[i] - 1))
      })),
      estimate = eta_B
    )
  }
  
  # stopifnot_(
  #   checkmate::test_numeric(
  #     x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
  #   ),
  #   "Argument {.arg probs} must be a {.cls numeric} vector with values
  #     between 0 and 1."
  # )
  # p_i <- length(eta_pi)
  # p_s <- length(eta_A)
  # p_o <- length(eta_B)
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
  #   gamma_pi[paste0("q", 100 * probs[i])] <- eta_pi + q * sds[seq_len(p_i)]
  #   gamma_A[paste0("q", 100 * probs[i])] <- eta_A + q * sds[p_i + seq_len(p_s)]
  #   gamma_B[paste0("q", 100 * probs[i])] <- eta_B + q * sds[p_i + p_s + seq_len(p_o)]
  # }
  
  list(
    gamma_initial = gamma_pi, 
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
  eta_pi <- unlist(object$coefficients$eta_pi)
  K_i <- length(object$coef_names_initial)
  object$state_names <- unname(object$state_names)
  gamma_pi <- data.frame(
    state = unlist(lapply(object$state_names, function(x) x[-1])),
    parameter = rep(object$coef_names_initial, each = (S - 1)),
    estimate = unlist(eta_pi),
    cluster = rep(object$cluster_names, each = (S - 1) * K_i)
  )
  eta_A <- unlist(object$coefficients$eta_A)
  K_s <- length(object$coef_names_transition)
  gamma_A <- data.frame(
    state_from = unlist(object$state_names),
    state_to = rep(
      unlist(lapply(object$state_names, function(x) x[-1])), 
      each = S
    ),
    parameter = rep(object$coef_names_transition, each = S * (S - 1)),
    estimate = unlist(eta_A),
    cluster = rep(object$cluster_names, each = (S - 1) * S * K_s)
  )
  K_o <- length(object$coef_names_emission)
  if (object$n_channels == 1) {
    gamma_B <- data.frame(
      state = unlist(object$state_names),
      observations = rep(object$symbol_names[-1], each = S),
      parameter = rep(object$coef_names_emission, each = S * (M - 1)),
      estimate = unlist(object$coefficients$eta_B),
      cluster =  rep(object$cluster_names, each = S * (M - 1) * K_o)
    )
  } else {
    gamma_B <- data.frame(
      state = unlist(object$state_names),
      observations = rep(unlist(lapply(object$symbol_names, "[", -1)), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * (M[i] - 1))
      })),
      estimate = unlist(object$coefficients$eta_B),
      cluster =  unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$cluster_names, each = S * (M[i] - 1) * K_o)
      }))
    )
  }
  eta_omega <- c(object$coefficients$eta_omega)
  gamma_omega <- data.frame(
    cluster = object$cluster_names[-1],
    parameter = rep(object$coef_names_cluster, each = D - 1),
    estimate = eta_omega
  )
  
  # stopifnot_(
  #   checkmate::test_numeric(
  #     x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
  #   ),
  #   "Argument {.arg probs} must be a {.cls numeric} vector with values
  #     between 0 and 1."
  # )
  # p_i <- length(eta_pi)
  # p_s <- length(eta_A)
  # p_o <- length(eta_B)
  # p_d <- length(eta_omega)
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
  #   gamma_pi[paste0("q", 100 * probs[i])] <- eta_pi + q * sds[seq_len(p_i)]
  #   gamma_A[paste0("q", 100 * probs[i])] <- eta_A + q * sds[p_i + seq_len(p_s)]
  #   gamma_B[paste0("q", 100 * probs[i])] <- eta_B + q * sds[p_i + p_s + seq_len(p_o)]
  #   gamma_omega[paste0("q", 100 * probs[i])] <- eta_omega + q * sds[p_i + p_s + p_o + seq_len(p_d)]
  # }
  
  list(
    gamma_initial = gamma_pi, 
    gamma_transition = gamma_A, 
    gamma_emission = gamma_B,
    gamma_cluster = gamma_omega
  )
}
