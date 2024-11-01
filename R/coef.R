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
  coef_names <- attr(object$X_pi, "coef_names")
  sd_pi_X <- rep(
    c(
      if(coef_names[1] == "(Intercept)") 1 else NULL, 
      attr(object$X_pi, "X_sd")
    ), each = S
  )
  gamma_pi <- data.frame(
    state = object$state_names,
    parameter = rep(coef_names, each = S),
    estimate = c(object$gammas$pi) / sd_pi_X
  )
  coef_names <- attr(object$X_A, "coef_names")
  sd_A_X <- rep(
    c(
      if(coef_names[1] == "(Intercept)") 1 else NULL, 
      attr(object$X_A, "X_sd")
    ), each = S^2
  )
  gamma_A <- data.frame(
    state_from = object$state_names,
    state_to = rep(object$state_names, each = S),
    parameter = rep(coef_names, each = S^2),
    estimate = c(object$gammas$A) / sd_A_X
  )
  coef_names <- attr(object$X_B, "coef_names")
  sd_B_X <- c(
    if(coef_names[1] == "(Intercept)") 1 else NULL, 
    attr(object$X_B, "X_sd")
  )
  if (object$n_channels == 1) {
    gamma_B <- data.frame(
      state = object$state_names,
      observation = rep(object$symbol_names, each = S),
      parameter = rep(coef_names, each = S * M),
      estimate = c(object$gammas$B) / rep(sd_B_X, each = S * M)
    )
  } else {
    gamma_B <- do.call(
      rbind,
      lapply(
        seq_len(object$n_channels), function(i) {
          data.frame(
            state = object$state_names,
            observation = rep(object$symbol_names[[i]], each = S),
            parameter = rep(coef_names, each = S * M[i]),
            estimate = c(object$gammas$B[[i]]) / rep(sd_B_X, each = S * M[i])
          )
        }
      )
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
    q_pi <- fast_quantiles(matrix(unlist(object$boot$gamma_pi), ncol = B), probs)
    q_A <- fast_quantiles(matrix(unlist(object$boot$gamma_A), ncol = B), probs)
    q_B <- fast_quantiles(matrix(unlist(object$boot$gamma_B), ncol = B), probs)
    sd_B_X <- unlist(lapply(seq_along(M),
                            function(i) rep(sd_B_X, each = S * M[i])))
    for(i in seq_along(probs)) {
      gamma_pi[paste0("q", 100 * probs[i])] <- q_pi[, i] / sd_pi_X
      gamma_A[paste0("q", 100 * probs[i])] <- q_A[, i] / sd_A_X
      gamma_B[paste0("q", 100 * probs[i])] <- q_B[, i] / sd_B_X
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
  object$state_names <- unname(object$state_names)
  coef_names <- attr(object$X_pi, "coef_names")
  sd_pi_X <- rep(
    c(
      if(coef_names[1] == "(Intercept)") 1 else NULL, 
      attr(object$X_pi, "X_sd")
    ), each = S
  )
  K_pi <- length(coef_names)
  gamma_pi <- data.frame(
    cluster = rep(object$cluster_names, each = S * K_pi),
    state = unlist(object$state_names),
    parameter = rep(coef_names, each = S),
    estimate = unlist(object$gammas$pi) / sd_pi_X
  )
  coef_names <- attr(object$X_A, "coef_names")
  sd_A_X <- rep(
    c(
      if(coef_names[1] == "(Intercept)") 1 else NULL, 
      attr(object$X_A, "X_sd")
    ), each = S^2
  )
  K_A <- length(coef_names)
  gamma_A <- data.frame(
    cluster = rep(object$cluster_names, each = S * S * K_A),
    state_from = unlist(object$state_names),
    state_to = rep(unlist(object$state_names), each = S),
    parameter = rep(coef_names, each = S * S),
    estimate = unlist(object$gammas$A) / sd_A_X
  )
  coef_names <- attr(object$X_B, "coef_names")
  sd_B_X <- c(
    if(coef_names[1] == "(Intercept)") 1 else NULL, 
    attr(object$X_B, "X_sd")
  )
  K_B <- length(coef_names)
  if (object$n_channels == 1) {
    gamma_B <- data.frame(
      cluster =  rep(object$cluster_names, each = S * M * K_B),
      state = unlist(object$state_names),
      observation = rep(object$symbol_names, each = S),
      parameter = rep(coef_names, each = S * M),
      estimate = unlist(object$gammas$B) / rep(sd_B_X, each = S * M)
    )
  } else {
    gamma_B <- do.call(
      rbind, 
      lapply(
        seq_len(object$n_clusters), function(d) {
          do.call(
            rbind,
            lapply(
              seq_len(object$n_channels), function(i) {
                data.frame(
                  cluster = rep(object$cluster_names[d], each = S * M[i] * K_B),
                  state = object$state_names[[d]],
                  observation = rep(object$symbol_names[[i]], each = S),
                  parameter = rep(coef_names, each = S * M[i]),
                  estimate = c(object$gammas$B[[d]][[i]]) / rep(sd_B_X, each = S * M[i])
                )
              }
            )
          )
        }
      )
    )
  }
  coef_names <- attr(object$X_omega, "coef_names")
  sd_omega_X <- rep(
    c(
      if(coef_names[1] == "(Intercept)") 1 else NULL, 
      attr(object$X_omega, "X_sd")
    ), each = D
  )
  gamma_omega <- data.frame(
    cluster = object$cluster_names,
    parameter = rep(coef_names, each = D),
    estimate = c(object$gammas$omega) / sd_omega_X
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
    q_pi <- fast_quantiles(matrix(unlist(object$boot$gamma_pi), ncol = B), probs)
    q_A <- fast_quantiles(matrix(unlist(object$boot$gamma_A), ncol = B), probs)
    q_B <- fast_quantiles(matrix(unlist(object$boot$gamma_B), ncol = B), probs)
    q_omega <- fast_quantiles(matrix(unlist(object$boot$gamma_omega), ncol = B), probs)
    sd_B_X <- unlist(lapply(seq_along(M),
                            function(i) rep(sd_B_X, each = S * M[i])))
    for(i in seq_along(probs)) {
      gamma_pi[paste0("q", 100 * probs[i])] <- q_pi[, i] / sd_pi_X
      gamma_A[paste0("q", 100 * probs[i])] <- q_A[, i] / sd_A_X
      gamma_B[paste0("q", 100 * probs[i])] <- q_B[, i] / sd_B_X
      gamma_omega[paste0("q", 100 * probs[i])] <- q_omega[, i] / sd_omega_X
    }
  }
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B,
    cluster = gamma_omega
  )
}
