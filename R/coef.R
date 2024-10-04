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
  gamma_pi <- data.frame(
    state = object$state_names,
    parameter = rep(object$coef_names_initial, each = S),
    estimate = c(object$gammas$pi)
  )
  gamma_A <- data.frame(
    state_from = object$state_names,
    state_to = rep(object$state_names, each = S),
    parameter = rep(object$coef_names_transition, each = S^2),
    estimate = c(object$gammas$A)
  )
  if (object$n_channels == 1) {
    gamma_B <- data.frame(
      state = object$state_names,
      observation = rep(object$symbol_names, each = S),
      parameter = rep(object$coef_names_emission, each = S * M),
      estimate = c(object$gammas$B)
    )
  } else {
    gamma_B <- data.frame(
      state = object$state_names,
      observation = rep(unlist(object$symbol_names), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * M[i])
      })),
      estimate = unlist(object$gammas$B)
    )
  }
  
  # stopifnot_(
  #   checkmate::test_numeric(
  #     x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
  #   ),
  #   "Argument {.arg probs} must be a {.cls numeric} vector with values
  #     between 0 and 1."
  # )
  # p_i <- length(gamma_pi)
  # p_s <- length(gamma_A)
  # p_o <- length(gamma_B)
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
  #   gamma_pi[paste0("q", 100 * probs[i])] <- gamma_pi + q * sds[seq_len(p_i)]
  #   gamma_A[paste0("q", 100 * probs[i])] <- gamma_A + q * sds[p_i + seq_len(p_s)]
  #   gamma_B[paste0("q", 100 * probs[i])] <- gamma_B + q * sds[p_i + p_s + seq_len(p_o)]
  # }
  
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B
  )
}
#' @rdname coef
#' @export
coef.mnhmm <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  
  S <- object$n_states
  M <- object$n_symbols
  D <- object$n_clusters
  K_i <- length(object$coef_names_initial)
  object$state_names <- unname(object$state_names)
  gamma_pi <- data.frame(
    state = unlist(object$state_names),
    parameter = rep(object$coef_names_initial, each = S),
    estimate = unlist(object$gammas$pi),
    cluster = rep(object$cluster_names, each = S * K_i)
  )
  K_s <- length(object$coef_names_transition)
  gamma_A <- data.frame(
    state_from = unlist(object$state_names),
    state_to = rep(unlist(object$state_names), each = S),
    parameter = rep(object$coef_names_transition, each = S * S),
    estimate = unlist(object$gammas$A),
    cluster = rep(object$cluster_names, each = S * S * K_s)
  )
  K_o <- length(object$coef_names_emission)
  if (object$n_channels == 1) {
    gamma_B <- data.frame(
      state = unlist(object$state_names),
      observations = rep(object$symbol_names, each = S),
      parameter = rep(object$coef_names_emission, each = S * M),
      estimate = unlist(object$gammas$B),
      cluster =  rep(object$cluster_names, each = S * M * K_o)
    )
  } else {
    gamma_B <- data.frame(
      state = unlist(object$state_names),
      observations = rep(unlist(object$symbol_names), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * M[i])
      })),
      estimate = unlist(object$gammas$B),
      cluster =  unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$cluster_names, each = S * M[i] * K_o)
      }))
    )
  }
  gamma_omega <- data.frame(
    cluster = object$cluster_names,
    parameter = rep(object$coef_names_cluster, each = D),
    estimate = c(object$gammas$omega)
  )
  
  # stopifnot_(
  #   checkmate::test_numeric(
  #     x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
  #   ),
  #   "Argument {.arg probs} must be a {.cls numeric} vector with values
  #     between 0 and 1."
  # )
  # p_i <- length(gamma_pi)
  # p_s <- length(gamma_A)
  # p_o <- length(gamma_B)
  # p_d <- length(gamma_omega)
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
  #   gamma_pi[paste0("q", 100 * probs[i])] <- gamma_pi + q * sds[seq_len(p_i)]
  #   gamma_A[paste0("q", 100 * probs[i])] <- gamma_A + q * sds[p_i + seq_len(p_s)]
  #   gamma_B[paste0("q", 100 * probs[i])] <- gamma_B + q * sds[p_i + p_s + seq_len(p_o)]
  #   gamma_omega[paste0("q", 100 * probs[i])] <- gamma_omega + q * sds[p_i + p_s + p_o + seq_len(p_d)]
  # }
  
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B,
    cluster = gamma_omega
  )
}
