#' Get the Estimated Initial, Transition, and Emission Probabilities for NHMM 
#' or MNHMM
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param nsim Non-negative integer defining the number of samples from the 
#' normal approximation of the model parameters used in 
#' computing the approximate quantiles of the estimates. If `0`, only point 
#' estimates are returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.975)`.
#' @rdname get_probs
#' @export
get_probs <- function(model, ...) {
  UseMethod("get_probs", model)
}
#' @rdname get_probs
#' @export
get_probs.nhmm <- function(model, data = NULL, nsim = 0, 
                           probs = c(0.025, 0.975)) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  if (!is.null(data)) {
    model <- update(model, data = data)
  }
  S <- model$n_states
  M <- model$n_symbols
  C <- model$n_channels
  beta_i_raw <- stan_to_cpp_initial(
    model$estimation_results$parameters$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    model$estimation_results$parameters$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    model$estimation_results$parameters$beta_o_raw,
    1,
    C > 1
  )
  if (intercept_only(model$initial_formula)) {
    X_initial <- t(model$X_initial[1, , drop = FALSE])
    ids_pi <- "all"
  } else {
    X_initial <- t(model$X_initial)
    ids_pi <- seq_len(model$n_sequences)
  }
  if (intercept_only(model$transition_formula)) {
    X_transition <- aperm(model$X_transition[1, 1, , drop = FALSE], c(3, 1, 2))
    ids_A <- "all"
    times_A <- "all"
  } else {
    X_transition <- aperm(model$X_transition, c(3, 1, 2))
    ids_A <- seq_len(model$n_sequences)
    times_A <- colnames(model$observations)
  }
  if (intercept_only(model$emission_formula)) {
    X_emission <- aperm(model$X_emission[1, 1, , drop = FALSE], c(3, 1, 2))
    ids_B <- "all"
    times_B <- "all"
  } else {
    X_emission <- aperm(model$X_emission, c(3, 1, 2))
    ids_B <- seq_len(model$n_sequences)
    times_B <- colnames(model$observations)
  }
  initial_probs <- get_pi(beta_i_raw, X_initial, 0)
  transition_probs <- get_A(beta_s_raw, X_transition, 0)
  emission_probs <- if (model$n_channels == 1) {
    get_B(beta_o_raw, X_emission, 0) 
  } else {
    get_multichannel_B(beta_o_raw, X_emission, S, C, M, 0) 
  }
  initial_probs <- data.frame(
    id = rep(ids_pi, each = S),
    state = model$state_names,
    estimate = c(initial_probs)
  )
  transition_probs <- data.frame(
    id = rep(ids_A, each = S^2),
    time = rep(times_A, each = S^2 * length(ids_A)),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = unlist(transition_probs)
  )
  emission_probs <- data.frame(
    id = rep(ids_B, each = S * M),
    time = rep(times_B, each = S * M * length(ids_B)),
    state = model$state_names,
    observation = rep(model$symbol_names, each = S),
    estimate = unlist(emission_probs)
  )
  if (nsim > 0) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    chol_precision <- chol(-model$estimation$hessian)
    U <- backsolve(chol_precision, diag(ncol(chol_precision)))
    x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
    x <- t(sweep(x, 2, c(beta_i_raw, beta_s_raw, beta_o_raw), "+"))
    p_i <- length(beta_i_raw)
    p_s <- length(beta_s_raw)
    p_o <- length(beta_o_raw)
    
    samples <- apply(
      x[seq_len(p_i), ], 2, function(z) {
        z <- array(z, dim = dim(beta_i_raw))
        get_pi(z, X_initial)
      }
    )
    quantiles <- fast_quantiles(samples, probs)
    for(i in seq_along(probs)) {
      initial_probs[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    samples <- apply(
      x[p_i + seq_len(p_s), ], 2, function(z) {
        z <- array(z, dim = dim(beta_s_raw))
        unlist(get_A(aperm(z, c(2, 3, 1)), X_transition))
      }
    )
    quantiles <- fast_quantiles(samples, probs)
    for(i in seq_along(probs)) {
      transition_probs[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    samples <- apply(
      x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
        z <- array(z, dim = dim(beta_o_raw))
        unlist(get_B(aperm(z, c(2, 3, 1)), X_emission))
      }
    )
    quantiles <- fast_quantiles(samples, probs)
    for(i in seq_along(probs)) {
      emission_probs[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
  }
  list(
    initial_probs = initial_probs, 
    transition_probs = 
      remove_voids(transition_probs, time, id, transition_probs),
    emission_probs = remove_voids(emission_probs, time, id, emission_probs)
  )
}
#' @rdname get_probs
#' @export
get_probs.mnhmm <- function(model, data = NULL, nsim = 0, 
                            probs = c(0.025, 0.975)) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  
  if (!is.null(data)) {
    model <- update(model, data = data)
  }
  S <- model$n_states
  M <- model$n_symbols
  C <- model$n_channels
  beta_i_raw <- stan_to_cpp_initial(
    model$estimation_results$parameters$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    model$estimation_results$parameters$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    model$estimation_results$parameters$beta_o_raw,
    1,
    C > 1
  )
  theta_raw <- model$estimation_results$parameters$theta_raw
  
  if (intercept_only(model$initial_formula)) {
    X_initial <- t(model$X_initial[1, , drop = FALSE])
    ids_pi <- "all"
  } else {
    X_initial <- t(model$X_initial)
    ids_pi <- seq_len(model$n_sequences)
  }
  if (intercept_only(model$transition_formula)) {
    X_transition <- aperm(model$X_transition[1, 1, , drop = FALSE], c(3, 1, 2))
    ids_A <- "all"
    times_A <- "all"
  } else {
    X_transition <- aperm(model$X_transition, c(3, 1, 2))
    ids_A <- seq_len(model$n_sequences)
    times_A <- colnames(model$observations)
  }
  if (intercept_only(model$emission_formula)) {
    X_emission <- aperm(model$X_emission[1, 1, , drop = FALSE], c(3, 1, 2))
    ids_B <- "all"
    times_B <- "all"
  } else {
    X_emission <- aperm(model$X_emission, c(3, 1, 2))
    ids_B <- seq_len(model$n_sequences)
    times_B <- colnames(model$observations)
  }
  if (intercept_only(model$cluster_formula)) {
    X_cluster <- t(model$X_cluster[1, , drop = FALSE])
    ids_theta <- "all"
  } else {
    X_cluster <- t(model$X_cluster)
    ids_theta <- seq_len(model$n_sequences)
  }
  initial_probs <- get_pi(beta_i_raw, X_initial, 0)
  transition_probs <- get_A(beta_s_raw, X_transition, 0)
  emission_probs <- if (model$n_channels == 1) {
    get_B(beta_o_raw, X_emission, 0) 
  } else {
    get_multichannel_B(beta_o_raw, X_emission, S, C, M, 0) 
  }
  cluster_probs <- get_theta(theta_raw, X_cluster, 0)
  initial_probs <- data.frame(
    id = rep(ids_pi, each = S),
    state = model$state_names,
    estimate = c(initial_probs)
  )
  transition_probs <- data.frame(
    id = rep(ids_A, each = S^2),
    time = rep(times_A, each = S^2 * length(ids_A)),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = unlist(transition_probs)
  )
  emission_probs <- data.frame(
    id = rep(ids_B, each = S * M),
    time = rep(times_B, each = S * M * length(ids_B)),
    state = model$state_names,
    observation = rep(model$symbol_names, each = S),
    estimate = unlist(emission_probs)
  )
  if (nsim > 0) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    chol_precision <- chol(-model$estimation$hessian)
    U <- backsolve(chol_precision, diag(ncol(chol_precision)))
    x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
    x <- t(sweep(x, 2, c(beta_i_raw, beta_s_raw, beta_o_raw), "+"))
    p_i <- length(beta_i_raw)
    p_s <- length(beta_s_raw)
    p_o <- length(beta_o_raw)
    
    samples <- apply(
      x[seq_len(p_i), ], 2, function(z) {
        z <- array(z, dim = dim(beta_i_raw))
        get_pi(z, X_initial)
      }
    )
    quantiles <- fast_quantiles(samples, probs)
    for(i in seq_along(probs)) {
      initial_probs[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    samples <- apply(
      x[p_i + seq_len(p_s), ], 2, function(z) {
        z <- array(z, dim = dim(beta_s_raw))
        unlist(get_A(aperm(z, c(2, 3, 1)), X_transition))
      }
    )
    quantiles <- fast_quantiles(samples, probs)
    for(i in seq_along(probs)) {
      transition_probs[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
    samples <- apply(
      x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
        z <- array(z, dim = dim(beta_o_raw))
        unlist(get_B(aperm(z, c(2, 3, 1)), X_emission))
      }
    )
    quantiles <- fast_quantiles(samples, probs)
    for(i in seq_along(probs)) {
      emission_probs[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
  }
  list(
    initial_probs = initial_probs, 
    transition_probs = 
      remove_voids(transition_probs, time, id, transition_probs),
    emission_probs = remove_voids(emission_probs, time, id, emission_probs),
    cluster_probs = cluster_probs
  )
}