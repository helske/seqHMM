#' Get the Estimated Regression Coefficients of Non-Homogeneous Hidden Markov 
#' Models
#' 
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param probs Vector defining the quantiles of interest. When `NULL` (default), 
#' no quantiles are computed. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param ... Ignored.
#' @rdname coef
#' @export
coef.nhmm <- function(object, probs = NULL, ...) {
  S <- object$n_states
  M <- object$n_symbols
  coef_names_pi <- attr(object$X_pi, "coef_names")
  R_inv_pi <- attr(object$X_pi, "R_inv")
  X_mean_pi <- attr(object$X_pi, "X_mean")
  gamma_pi <- gamma_std_to_gamma(
    object$gammas$pi, R_inv_pi, coef_names_pi, X_mean_pi
  )
  gamma_pi <- data.table(
    state = object$state_names,
    parameter = rep(coef_names_pi, each = S),
    estimate = c(gamma_pi)
  )
  coef_names_A <- attr(object$X_A, "coef_names")
  R_inv_A <- attr(object$X_A, "R_inv")
  X_mean_A <- attr(object$X_A, "X_mean")
  K <- length(coef_names_A)
  gamma_A <- gamma_std_to_gamma(
    object$gammas$A, R_inv_A, coef_names_A, X_mean_A
  )
  gamma_A <- data.table(
    state_from = rep(object$state_names, each = S * K),
    state_to = object$state_names,
    parameter = rep(coef_names_A, each = S),
    estimate = c(gamma_A)
  )
  coef_names_B <- attr(object$X_B, "coef_names")
  R_inv_B <- attr(object$X_B, "R_inv")
  X_mean_B <- attr(object$X_B, "X_mean")
  K <- length(coef_names_B)
  gamma_B <- gamma_std_to_gamma(
    object$gammas$B, R_inv_B, coef_names_B, X_mean_B
  )
  if (object$n_channels == 1) {
    gamma_B <- data.table(
      state = rep(object$state_names, each = M * K),
      observation = object$symbol_names,
      parameter = rep(coef_names_B, each = M),
      estimate = c(gamma_B)
    )
  } else {
    gamma_B <- do.call(
      rbind,
      lapply(
        seq_len(object$n_channels), function(i) {
          data.table(
            state = rep(object$state_names, each = M[i] * K),
            observation = object$symbol_names[[i]],
            parameter = rep(coef_names_B, each = M[i]),
            estimate = c(gamma_B[[i]])
          )
        }
      )
    )
  }
  if (!is.null(probs)) {
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
    nsim <- length(object$boot$gamma_pi)
    boot_gamma_pi <- lapply(
      object$boot$gamma_pi, 
      gamma_std_to_gamma, R_inv = R_inv_pi,
      coef_names = coef_names_pi, X_mean = X_mean_pi
    )
    boot_gamma_A <- lapply(
      object$boot$gamma_A, 
      gamma_std_to_gamma, R_inv = R_inv_A,
      coef_names = coef_names_A, X_mean = X_mean_A
    )
    boot_gamma_B <- lapply(
      object$boot$gamma_B, 
      gamma_std_to_gamma, R_inv = R_inv_B,
      coef_names = coef_names_B, X_mean = X_mean_B
    )
    q_pi <- fast_quantiles(matrix(unlist(boot_gamma_pi), ncol = nsim), probs)
    q_A <- fast_quantiles(matrix(unlist(boot_gamma_A), ncol = nsim), probs)
    q_B <- fast_quantiles(matrix(unlist(boot_gamma_B), ncol = nsim), probs)
    probs <- paste0("q", 100 * probs)
    gamma_pi[, (probs) := data.table(q_pi)]
    gamma_A[, (probs) := data.table(q_A)]
    gamma_B[, (probs) := data.table(q_B)]
  }
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B
  )
}
#' @rdname coef
#' @export
coef.mnhmm <- function(object, probs = NULL, ...) {
  
  S <- object$n_states
  M <- object$n_symbols
  D <- object$n_clusters
  object$state_names <- unname(object$state_names)
  coef_names_pi <- attr(object$X_pi, "coef_names")
  K <- length(coef_names_pi)
  R_inv_pi <- attr(object$X_pi, "R_inv")
  X_mean_pi <- attr(object$X_pi, "X_mean")
  gamma_pi <- gamma_std_to_gamma(
    object$gammas$pi, R_inv_pi, coef_names_pi, X_mean_pi
  )
  gamma_pi <- data.table(
    cluster = rep(object$cluster_names, each = S * K),
    state = unlist(object$state_names),
    parameter = rep(coef_names_pi, each = S),
    estimate = unlist(gamma_pi)
  )
  coef_names_A <- attr(object$X_A, "coef_names")
  R_inv_A <- attr(object$X_A, "R_inv")
  X_mean_A <- attr(object$X_A, "X_mean")
  K <- length(coef_names_A)
  gamma_A <- gamma_std_to_gamma(
    object$gammas$A, R_inv_A, coef_names_A, X_mean_A
  )
  gamma_A <- data.table(
    cluster = rep(object$cluster_names, each = S * S * K),
    state_from = rep(unlist(object$state_names), each = S * K),
    state_to = unlist(object$state_names),
    parameter = rep(coef_names_A, each = S),
    estimate = unlist(gamma_A)
  )
  
  coef_names_B <- attr(object$X_B, "coef_names")
  R_inv_B <- attr(object$X_B, "R_inv")
  X_mean_B <- attr(object$X_B, "X_mean")
  K <- length(coef_names_B)
  gamma_B <- gamma_std_to_gamma(
    object$gammas$B, R_inv_B, coef_names_B, X_mean_B
  )
  if (object$n_channels == 1) {
    gamma_B <- data.table(
      cluster =  rep(object$cluster_names, each = S * M * K),
      state = rep(unlist(object$state_names), each = M * K),
      observation = object$symbol_names,
      parameter = rep(coef_names_B, each = M),
      estimate = unlist(gamma_B)
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
                data.table(
                  cluster = rep(object$cluster_names[d], each = S * M[i] * K),
                  state = rep(object$state_names[[d]], each = M[i] * K),
                  observation = object$symbol_names[[i]],
                  parameter = rep(coef_names_B, each = M[i]),
                  estimate = c(gamma_B[[d]][[i]])
                )
              }
            )
          )
        }
      )
    )
  }
  coef_names_omega <- attr(object$X_omega, "coef_names")
  R_inv_omega <- attr(object$X_omega, "R_inv")
  X_mean_omega <- attr(object$X_omega, "X_mean")
  gamma_omega <- gamma_std_to_gamma(
    object$gammas$omega, R_inv_omega, coef_names_omega, X_mean_omega
  )
  gamma_omega <- data.table(
    cluster = object$cluster_names,
    parameter = rep(coef_names_omega, each = D),
    estimate = c(gamma_omega)
  )
  if (!is.null(probs)) {
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
    nsim <- length(object$boot$gamma_pi)
    boot_gamma_pi <- lapply(
      object$boot$gamma_pi, 
      gamma_std_to_gamma, R_inv = R_inv_pi,
      coef_names = coef_names_pi, X_mean = X_mean_pi
    )
    boot_gamma_A <- lapply(
      object$boot$gamma_A, 
      gamma_std_to_gamma, R_inv = R_inv_A,
      coef_names = coef_names_A, X_mean = X_mean_A
    )
    boot_gamma_B <- lapply(
      object$boot$gamma_B, 
      gamma_std_to_gamma, R_inv = R_inv_B,
      coef_names = coef_names_B, X_mean = X_mean_B
    )
    boot_gamma_omega <- lapply(
      object$boot$gamma_omega, 
      gamma_std_to_gamma, R_inv = R_inv_omega,
      coef_names = coef_names_omega, X_mean = X_mean_omega
    )
    q_pi <- fast_quantiles(matrix(unlist(boot_gamma_pi), ncol = nsim), probs)
    q_A <- fast_quantiles(matrix(unlist(boot_gamma_A), ncol = nsim), probs)
    q_B <- fast_quantiles(matrix(unlist(boot_gamma_B), ncol = nsim), probs)
    q_omega <- fast_quantiles(matrix(unlist(boot_gamma_omega), ncol = nsim), probs)
    
    probs <- paste0("q", 100 * probs)
    gamma_pi[, (probs) := data.table(q_pi)]
    gamma_A[, (probs) := data.table(q_A)]
    gamma_B[, (probs) := data.table(q_B)]
    gamma_omega[, (probs) := data.table(q_omega)]
  }
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B,
    cluster = gamma_omega
  )
}
