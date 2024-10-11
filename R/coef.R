#' Get the Estimated Regression Coefficients of Non-Homogeneous Hidden Markov 
#' Models
#' 
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param probs Vector defining the quantiles of interest. When missing (default), 
#' no quantiles are computed. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param ... Ignored.
#' @rdname coef
#' @export
coef.nhmm <- function(object, probs, ...) {
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
  if (!missing(probs)) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    p_i <- length(unlist(object$gammas$pi))
    p_s <- length(unlist(object$gammas$A))
    p_o <- length(unlist(object$gammas$B))
    stopifnot_(
      !is.null(object$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    B <- length(object$boot$gamma_pi)
    q_pi <- seqHMM:::fast_quantiles(matrix(unlist(object$boot$gamma_pi), ncol = B), probs)
    q_A <- seqHMM:::fast_quantiles(matrix(unlist(object$boot$gamma_A), ncol = B), probs)
    q_B <- seqHMM:::fast_quantiles(matrix(unlist(object$boot$gamma_B), ncol = B), probs)
    for(i in seq_along(probs)) {
      gamma_pi[paste0("q", 100 * probs[i])] <- q_pi[, i]
      gamma_A[paste0("q", 100 * probs[i])] <- q_A[, i]
      gamma_B[paste0("q", 100 * probs[i])] <- q_B[, i]
    }
  }
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B
  )
}
#' @rdname coef
#' @export
coef.mnhmm <- function(object, probs, ...) {
  
  S <- object$n_states
  M <- object$n_symbols
  D <- object$n_clusters
  K_i <- length(object$coef_names_initial)
  object$state_names <- unname(object$state_names)
  gamma_pi <- data.frame(
    cluster = rep(object$cluster_names, each = S * K_i),
    state = unlist(object$state_names),
    parameter = rep(object$coef_names_initial, each = S),
    estimate = unlist(object$gammas$pi)
  )
  K_s <- length(object$coef_names_transition)
  gamma_A <- data.frame(
    cluster = rep(object$cluster_names, each = S * S * K_s),
    state_from = unlist(object$state_names),
    state_to = rep(unlist(object$state_names), each = S),
    parameter = rep(object$coef_names_transition, each = S * S),
    estimate = unlist(object$gammas$A)
  )
  K_o <- length(object$coef_names_emission)
  if (object$n_channels == 1) {
    gamma_B <- data.frame(
      cluster =  rep(object$cluster_names, each = S * M * K_o),
      state = unlist(object$state_names),
      observations = rep(object$symbol_names, each = S),
      parameter = rep(object$coef_names_emission, each = S * M),
      estimate = unlist(object$gammas$B)
    )
  } else {
    gamma_B <- data.frame(
      cluster =  unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$cluster_names, each = S * M[i] * K_o)
      })),
      state = unlist(object$state_names),
      observations = rep(unlist(object$symbol_names), each = S),
      parameter = unlist(lapply(seq_len(object$n_channels), function(i) {
        rep(object$coef_names_emission, each = S * M[i])
      })),
      estimate = unlist(object$gammas$B)
    )
  }
  gamma_omega <- data.frame(
    cluster = object$cluster_names,
    parameter = rep(object$coef_names_cluster, each = D),
    estimate = c(object$gammas$omega)
  )
  if(!missing(probs)) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    p_i <- length(unlist(object$gammas$pi))
    p_s <- length(unlist(object$gammas$A))
    p_o <- length(unlist(object$gammas$B))
    p_d <- length(object$gammas$omega)
    
    stopifnot_(
      !is.null(object$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    B <- length(object$boot$gamma_pi)
    q_pi <- seqHMM:::fast_quantiles(matrix(unlist(object$boot$gamma_pi), ncol = B), probs)
    q_A <- seqHMM:::fast_quantiles(matrix(unlist(object$boot$gamma_A), ncol = B), probs)
    q_B <- seqHMM:::fast_quantiles(matrix(unlist(object$boot$gamma_B), ncol = B), probs)
    q_omega <- seqHMM:::fast_quantiles(matrix(unlist(object$boot$gamma_omega), ncol = B), probs)
    for(i in seq_along(probs)) {
      gamma_pi[paste0("q", 100 * probs[i])] <- q_pi[, i]
      gamma_A[paste0("q", 100 * probs[i])] <- q_A[, i]
      gamma_B[paste0("q", 100 * probs[i])] <- q_B[, i]
      gamma_omega[paste0("q", 100 * probs[i])] <- q_omega[, i]
    }
  }
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B,
    cluster = gamma_omega
  )
}
