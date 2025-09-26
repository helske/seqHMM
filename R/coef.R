#' Get the Estimated Regression Coefficients of Non-Homogeneous Hidden Markov 
#' Models
#' 
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param probs Vector defining the quantiles of interest. When `NULL` (default), 
#' no quantiles are computed. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param ... Ignored.
#' @rdname coef
#' @return A list of data tables with the estimated coefficients for initial,
#' transition, emission (separate `data.table` for each response), and cluster 
#' probabilities (in case of mixture model). 
#' @export
coef.nhmm <- function(object, probs = NULL, ...) {
  S <- object$n_states
  M <- object$n_symbols
  states <- object$state_names
  symbols <- object$symbol_names
  responses <- object$responses
  coef_names_pi <- attr(object$X_pi, "coef_names")
  R_inv_pi <- attr(object$X_pi, "R_inv")
  X_mean_pi <- attr(object$X_pi, "X_mean")
  gamma_pi <- gamma_std_to_gamma(
    object$gammas$gamma_pi, R_inv_pi, coef_names_pi, X_mean_pi
  )
  gamma_pi <- data.table(
    state = states,
    coefficient = rep(coef_names_pi, each = S),
    estimate = c(gamma_pi)
  )
  coef_names_A <- attr(object$X_A, "coef_names")
  R_inv_A <- attr(object$X_A, "R_inv")
  X_mean_A <- attr(object$X_A, "X_mean")
  K <- length(coef_names_A)
  gamma_A <- gamma_std_to_gamma(
    object$gammas$gamma_A, R_inv_A, coef_names_A, X_mean_A
  )
  gamma_A <- data.table(
    state_from = rep(states, each = S * K),
    state_to = states,
    coefficient = rep(coef_names_A, each = S),
    estimate = c(gamma_A)
  )
  coef_names_B <- lapply(responses, \(y) attr(object$X_B[[y]], "coef_names"))
  R_inv_B <- lapply(responses, \(y) attr(object$X_B[[y]], "R_inv"))
  X_mean_B <- lapply(responses, \(y) attr(object$X_B[[y]], "X_mean"))
  gamma_B <- gamma_std_to_gamma(
    object$gammas$gamma_B, R_inv_B, coef_names_B, X_mean_B
  )
  .funB <- function(gammas, states, symbols, coef_names) {
    K <- length(coef_names)
    M <- length(symbols)
    data.table(
      state = rep(states, each = M * K),
      observation = symbols,
      coefficient = rep(coef_names, each = M),
      estimate = c(gammas)
    )
  }
  gamma_B <- stats::setNames(
    lapply(
      seq_along(responses), \(i) .funB(
        gamma_B[[i]], states, symbols[[i]], coef_names_B[[i]]
      )
    ),
    responses
  )
  
  if (!is.null(probs)) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    stopifnot_(
      !is.null(object$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    nsim <- length(object$boot$gamma_pi)
    p <- paste0("q", 100 * probs)
    boot_gamma_pi <- lapply(
      seq_along(object$boot$gamma_pi), 
      \(i) gamma_std_to_gamma(
        object$boot$gamma_pi[[i]], R_inv_pi, coef_names_pi, X_mean_pi
      )
    )
    q_pi <- fast_quantiles(matrix(unlist(boot_gamma_pi), ncol = nsim), probs)
    gamma_pi[, (p) := data.table(q_pi)]
    boot_gamma_A <- lapply(
      seq_along(object$boot$gamma_A), 
      \(i) gamma_std_to_gamma(
        object$boot$gamma_A[[i]], R_inv_A, coef_names_A, X_mean_A
      )
    )
    q_A <- fast_quantiles(matrix(unlist(boot_gamma_A), ncol = nsim), probs)
    gamma_A[, (p) := data.table(q_A)]
    for (i in seq_along(responses)) {
      boot_gamma_B <- lapply(
        lapply(object$boot$gamma_B, "[[", i), 
        \(gamma) gamma_std_to_gamma(
          gamma, R_inv_B[[i]], coef_names_B[[i]], X_mean_B[[i]]
        )
      )
      q_B <- fast_quantiles(matrix(unlist(boot_gamma_B), ncol = nsim), probs)
      gamma_B[[responses[i]]][, (p) := data.table(q_B)]
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
coef.mnhmm <- function(object, probs = NULL, ...) {
  
  S <- object$n_states
  M <- object$n_symbols
  D <- object$n_clusters
  states <- object$state_names
  symbols <- object$symbol_names
  clusters <- object$cluster_names
  responses <- object$responses
  coef_names_pi <- attr(object$X_pi, "coef_names")
  K <- length(coef_names_pi)
  R_inv_pi <- attr(object$X_pi, "R_inv")
  X_mean_pi <- attr(object$X_pi, "X_mean")
  gamma_pi <- gamma_std_to_gamma(
    object$gammas$gamma_pi, R_inv_pi, coef_names_pi, X_mean_pi
  )
  gamma_pi <- data.table(
    cluster = rep(clusters, each = S * K),
    state = unlist(states),
    coefficient = rep(coef_names_pi, each = S),
    estimate = unlist(gamma_pi)
  )
  coef_names_A <- attr(object$X_A, "coef_names")
  R_inv_A <- attr(object$X_A, "R_inv")
  X_mean_A <- attr(object$X_A, "X_mean")
  K <- length(coef_names_A)
  gamma_A <- gamma_std_to_gamma(
    object$gammas$gamma_A, R_inv_A, coef_names_A, X_mean_A
  )
  gamma_A <- data.table(
    cluster = rep(clusters, each = S * S * K),
    state_from = rep(unlist(states), each = S * K),
    state_to = unlist(states),
    coefficient = rep(coef_names_A, each = S),
    estimate = unlist(gamma_A)
  )
  
  coef_names_B <- lapply(responses, \(y) attr(object$X_B[[y]], "coef_names"))
  R_inv_B <- lapply(responses, \(y) attr(object$X_B[[y]], "R_inv"))
  X_mean_B <- lapply(responses, \(y) attr(object$X_B[[y]], "X_mean"))

  .funB <- function(gammas, i, clusters, states, symbols, coef_names) {
    K <- length(coef_names[[i]])
    M <- length(symbols[[i]])
    do.call(
      rbind,
      lapply(
        seq_along(clusters), \(d) {
          data.table(
            cluster = clusters[d],
            state = rep(states[[d]], each = M * K),
            observation = symbols[[i]],
            coefficient = rep(coef_names[[i]], each = M),
            estimate = c(gammas[[d]][[i]])
          )
        }
      )
    )
  }
  gamma_B <- gamma_std_to_gamma(
    object$gammas$gamma_B, R_inv_B, coef_names_B, X_mean_B
  )
  gamma_B <- stats::setNames(
    lapply(
      seq_along(responses), \(i) .funB(
        gamma_B, i, clusters, states, symbols, coef_names_B
      )
    ),
    responses
  )
  
  coef_names_omega <- attr(object$X_omega, "coef_names")
  R_inv_omega <- attr(object$X_omega, "R_inv")
  X_mean_omega <- attr(object$X_omega, "X_mean")
  gamma_omega <- gamma_std_to_gamma(
    object$gammas$gamma_omega, R_inv_omega, coef_names_omega, X_mean_omega
  )
  gamma_omega <- data.table(
    cluster = object$cluster_names,
    coefficient = rep(coef_names_omega, each = D),
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
    stopifnot_(
      !is.null(object$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    nsim <- length(object$boot$gamma_pi)
    p <- paste0("q", 100 * probs)
    boot_gamma_pi <- lapply(
      object$boot$gamma_pi, 
      gamma_std_to_gamma, R_inv = R_inv_pi,
      coef_names = coef_names_pi, X_mean = X_mean_pi
    )
    q_pi <- fast_quantiles(matrix(unlist(boot_gamma_pi), ncol = nsim), probs)
    gamma_pi[, (p) := data.table(q_pi)]
    boot_gamma_A <- lapply(
      object$boot$gamma_A, 
      gamma_std_to_gamma, R_inv = R_inv_A,
      coef_names = coef_names_A, X_mean = X_mean_A
    )
    q_A <- fast_quantiles(matrix(unlist(boot_gamma_A), ncol = nsim), probs)
    gamma_A[, (p) := data.table(q_A)]
    boot_gamma_B <- lapply(
      object$boot$gamma_B, 
      \(gamma) gamma_std_to_gamma(
        gamma, R_inv_B, coef_names_B, X_mean_B
      )
    )
    for (i in seq_along(responses)) {
      gamma <- lapply(boot_gamma_B, \(x) {
        lapply(x, \(j) j[[i]])
      })
      q_B <- fast_quantiles(matrix(unlist(gamma), ncol = nsim), probs)
      gamma_B[[responses[i]]][, (p) := data.table(q_B)]
    }
    boot_gamma_omega <- lapply(
      object$boot$gamma_omega, 
      gamma_std_to_gamma, R_inv = R_inv_omega,
      coef_names = coef_names_omega, X_mean = X_mean_omega
    )
    q_omega <- fast_quantiles(matrix(unlist(boot_gamma_omega), ncol = nsim), probs)
    gamma_omega[, (p) := data.table(q_omega)]
  }
  list(
    initial = gamma_pi, 
    transition = gamma_A, 
    emission = gamma_B,
    cluster = gamma_omega
  )
}
