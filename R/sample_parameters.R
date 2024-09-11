#' Samples parameters of non-homogeneous Markov hidden models using the normal approximation.
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param nsim A non-negative integer defining the number of samples from the 
#' normal approximation.
#' @noRd
sample_parameters_nhmm <- function(model, nsim) {
  X_initial <- t(model$X_initial)
  X_transition <- aperm(model$X_transition, c(3, 1, 2))
  X_emission <- aperm(model$X_emission, c(3, 1, 2))
  
  chol_precision <- chol(-model$estimation$hessian)
  U <- backsolve(chol_precision, diag(ncol(chol_precision)))
  x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
  beta_i_raw <- model$coefficients$beta_i_raw
  beta_s_raw <- model$coefficients$beta_s_raw
  beta_o_raw <- model$coefficients$beta_o_raw
  pars <- c(beta_i_raw, beta_s_raw, beta_o_raw)
  p_i <- length(beta_i_raw)
  p_s <- length(beta_s_raw)
  p_o <- length(beta_o_raw)
  x <- t(sweep(x, 2, pars, "+"))
  samples_pi <- apply(
    x[seq_len(p_i), ], 2, function(z) {
      z <- array(z, dim = dim(beta_i_raw))
      get_pi(z, X_initial, 0)
    }
  )
  samples_A <- apply(
    x[p_i + seq_len(p_s), ], 2, function(z) {
      z <- array(z, dim = dim(beta_s_raw))
      unlist(get_A(stan_to_cpp_transition(z, D), X_transition, 0))
    }
  )
  if (model$n_channels == 1) {
    samples_B <- apply(
      x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
        z <- array(z, dim = dim(beta_o_raw))
        unlist(get_B(stan_to_cpp_emission(z, D, FALSE), X_emission, 0, 0))
      }
    )
  } else { 
    samples_B <- apply(
      x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
        z <- array(z, dim = dim(beta_o_raw))
        unlist(get_multichannel_B(
          stan_to_cpp_emission(z, D, TRUE), X_emission,
          model$n_states, model$n_channels, model$n_symbols, 0, 0
        ))
      }
    )
  }
  list(
    pi = samples_pi,
    A = samples_A,
    B = samples_B
  )
}

sample_parameters_mnhmm <- function(model, nsim) {
  X_initial <- t(model$X_initial)
  X_transition <- aperm(model$X_transition, c(3, 1, 2))
  X_emission <- aperm(model$X_emission, c(3, 1, 2))
  X_cluster <- t(model$X_cluster)
  
  chol_precision <- chol(-model$estimation$hessian)
  U <- backsolve(chol_precision, diag(ncol(chol_precision)))
  x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
  beta_i_raw <- model$coefficients$beta_i_raw
  beta_s_raw <- model$coefficients$beta_s_raw
  beta_o_raw <- model$coefficients$beta_o_raw 
  theta_raw <- model$coefficients$theta_raw
  pars <- c(beta_i_raw, beta_s_raw, beta_o_raw, theta_raw)
  D <- model$n_clusters
  p_i <- length(beta_i_raw)
  p_s <- length(beta_s_raw)
  p_o <- length(beta_o_raw)
  p_d <- length(theta_raw)
  x <- t(sweep(x, 2, pars, "+"))
  
  samples_pi <- apply(
    x[seq_len(p_i), ], 2, function(z) {
      b <- array(z, dim = dim(beta_i_raw))
      pi <- vector("list", D)
      for (d in seq_len(D)) {
        pi[[d]] <- get_pi(
          array(b[d, ,], dim = dim(b)[-1]), 
          X_initial, 0)
      }
      do.call("rbind", pi)
    }
  )
  samples_A <- apply(
    x[p_i + seq_len(p_s), ], 2, function(z) {
      b <- array(z, dim = dim(beta_s_raw))
      A <- vector("list", D)
      for (d in seq_len(D)) {
        A[[d]] <- unlist(get_A(
          stan_to_cpp_transition(array(b[d, , , ], dim = dim(b)[-1])), 
          X_transition, 0))
      }
      do.call("rbind", A)
    }
  )
  if (model$n_channels == 1) {
    samples_B <- apply(
      x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
        b <- array(z, dim = dim(beta_o_raw))
        B <- vector("list", D)
        for (d in seq_len(D)) {
          B[[d]] <- unlist(get_B(
            stan_to_cpp_emission(array(b[d, , , ], dim = dim(b)[-1]), 1, FALSE), 
            X_emission, 0, 0))
        }
        do.call("rbind", B)
      }
    )
  } else {
    samples_B <- apply(
      x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
        b <- array(z, dim = dim(beta_o_raw))
        B <- vector("list", D)
        for (d in seq_len(D)) {
          B[[d]] <- unlist(get_multichannel_B(
            stan_to_cpp_emission(array(b[d, ], dim = dim(b)[-1]), 1, TRUE), 
            X_emission, 0, 0))
        }
        do.call("rbind", B)
      }
    )
  }
  
  samples_omega <- apply(
    x[p_i + p_s + p_o + seq_len(p_d), ], 2, function(z) {
      z <- array(z, dim = dim(theta_raw))
      unlist(get_omega(z, X_cluster, 0))
    }
  )
  
  list(
    pi = samples_pi,
    A = samples_A,
    B = samples_B,
    omega = samples_omega
  )
}
#' Convert matrices of samples to data frames
#' @noRd
samples_to_df <- function(object, x, dontchange_colnames) {
  nsim <- ncol(x$pi)
  T_ <- object$length_of_sequences
  N <- object$n_sequences
  S <- object$n_states
  M <- object$n_symbols
  C <- object$n_channels
  D <- object$n_clusters
  mix <- D > 1
  if (C == 1) {
    ids <- rownames(object$observations)
    times <- colnames(object$observations)
    symbol_names <- list(object$symbol_names)
  } else {
    ids <- rownames(object$observations[[1]])
    times <- colnames(object$observations[[1]])
    symbol_names <- object$symbol_names
  }
  out <- vector("list", 3 + (D > 1))
  names(out) <- c("initial_probs", "transition_probs", "emission_probs", 
                  if (D > 1) "cluster_probs")
  out$initial_probs <- data.frame(
    cluster = if (mix) rep(object$cluster_names, each = S * N),
    id = rep(ids, each = S),
    state = object$state_names,
    estimate = c(x$pi),
    replication = rep(seq_len(nsim), each = nrow(x$pi))
  )
  out$transition_probs <- data.frame(
    if (mix) cluster = rep(object$cluster_names, each = S^2 * T_ * N),
    id = rep(ids, each = S^2 * T_),
    time = rep(times, each = S^2),
    state_from = object$state_names,
    state_to = rep(object$state_names, each = S),
    estimate = c(x$A),
    replication = rep(seq_len(nsim), each = nrow(x$A))
  )
  out$emission_probs <- data.frame(
    if (mix) cluster = rep(object$cluster_names, each = S * sum(M) * T_ * N),
    id = unlist(lapply(seq_len(C), function(i) rep(ids, each = S * M[i] * T_))),
    time = unlist(lapply(seq_len(C), function(i) rep(times, each = S * M[i]))),
    state = object$state_names,
    channel = rep(object$channel_names, S * M * T_ * N),
    observation = rep(unlist(symbol_names), each = S),
    estimate = c(x$B),
    replication = rep(seq_len(nsim), each = nrow(x$B))
  )
  if (mix) {
    names(out$initial_probs)[1] <- "cluster"
    names(out$transition_probs)[1] <- "cluster"
    names(out$emission_probs)[1] <- "cluster"
  }
  if (C == 1) out$emission_probs$channel <- NULL
  if (!dontchange_colnames) {
    colnames(out$initial_probs)[1 + mix] <- object$id_variable
    colnames(out$transition_probs)[1 + mix] <- object$id_variable
    colnames(out$transition_probs)[2 + mix] <- object$time_variable
    colnames(out$emission_probs)[1 + mix] <- object$id_variable
    colnames(out$emission_probs)[2 + mix] <- object$time_variable
  }
  if (mix) {
    out$cluster_probs <- data.frame(
      cluster = object$cluster_names,
      id = rep(ids, each = D),
      estimate = c(x$omega),
      replication = rep(seq_len(nsim), each = nrow(x$omega))
    )
    if (!dontchange_colnames) {
      colnames(out$cluster_probs)[2] <- object$id_variable
    }
  }
  out
}