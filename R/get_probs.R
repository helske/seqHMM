#' Get the Estimated Initial, Transition, and Emission Probabilities for NHMM 
#' or MNHMM
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param newdata An optional data frame containing the new data to be used in 
#' computing the probabilities.
#' @param nsim Non-negative integer defining the number of samples from the 
#' normal approximation of the model parameters used in 
#' computing the approximate quantiles of the estimates. If `0`, only point 
#' estimates are returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`.
#' @param ... Ignored.
#' @rdname get_probs
#' @export
get_probs <- function(model, ...) {
  UseMethod("get_probs", model)
}
#' @rdname get_probs
#' @export
get_probs.nhmm <- function(model, newdata = NULL, nsim = 0, 
                           probs = c(0.025, 0.5, 0.975), ...) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  if (!is.null(newdata)) {
    model <- update(model, newdata = newdata)
  }
  S <- model$n_states
  M <- model$n_symbols
  C <- model$n_channels
  N <- model$n_sequences
  T <- model$length_of_sequences
  beta_i_raw <- stan_to_cpp_initial(
    model$coefficients$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    model$coefficients$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    model$coefficients$beta_o_raw,
    1,
    C > 1
  )
  X_initial <- t(model$X_initial)
  X_transition <- aperm(model$X_transition, c(3, 1, 2))
  X_emission <- aperm(model$X_emission, c(3, 1, 2))
  initial_probs <- get_pi(beta_i_raw, X_initial, 0)
  transition_probs <- get_A(beta_s_raw, X_transition, 0)
  emission_probs <- if (C == 1) {
    get_B(beta_o_raw, X_emission, 0, 0) 
  } else {
    get_multichannel_B(beta_o_raw, X_emission, S, C, M, 0, 0) 
  } 
  if (C == 1) {
    ids <- rownames(model$observations)
    times <- colnames(model$observations)
    symbol_names <- list(model$symbol_names)
  } else {
    ids <- rownames(model$observations[[1]])
    times <- colnames(model$observations[[1]])
    symbol_names <- model$symbol_names
  }
  initial_probs <- data.frame(
    id = rep(ids, each = S),
    state = model$state_names,
    estimate = c(initial_probs)
  )
  colnames(initial_probs)[1] <- model$id_variable
  transition_probs <- data.frame(
    id = rep(ids, each = S^2 * T),
    time = rep(times, each = S^2),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = unlist(transition_probs)
  )
  colnames(transition_probs)[1] <- model$id_variable
  colnames(transition_probs)[2] <- model$time_variable
  emission_probs <- do.call(
    "rbind", 
    lapply(seq_len(C), function(i) {
      data.frame(
        id = rep(ids, each = S * M[i] * T),
        time = rep(times, each = S * M[i]),
        state = model$state_names,
        channel = model$channel_names[i],
        observation = rep(symbol_names[[i]], each = S),
        estimate = unlist(emission_probs[((i - 1) * N + 1):(i * N)])
      )
    })
  )
  colnames(emission_probs)[1] <- model$id_variable
  colnames(emission_probs)[2] <- model$time_variable
  
  if (nsim > 0) {
    out <- sample_parameters(model, nsim, probs)
    for(i in seq_along(probs)) {
      initial_probs[paste0("q", 100 * probs[i])] <- out$quantiles_pi[, i]
    }
    for(i in seq_along(probs)) {
      transition_probs[paste0("q", 100 * probs[i])] <- out$quantiles_A[, i]
    }
    for(i in seq_along(probs)) {
      emission_probs[paste0("q", 100 * probs[i])] <- out$quantiles_B[, i]
    }
  }
  rownames(initial_probs) <- NULL
  rownames(transition_probs) <- NULL
  rownames(emission_probs) <- NULL
  list(
    initial_probs = initial_probs, 
    transition_probs = remove_voids(model, transition_probs),
    emission_probs = remove_voids(model, emission_probs)
  )
}
#' @rdname get_probs
#' @export
get_probs.mnhmm <- function(model, newdata = NULL, nsim = 0, 
                            probs = c(0.025, 0.5, 0.975), ...) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  
  if (!is.null(newdata)) {
    model <- update(model, newdata = newdata)
  }
  
  T <- model$length_of_sequences
  N <- model$n_sequences
  S <- model$n_states
  M <- model$n_symbols
  C <- model$n_channels
  D <- model$n_clusters
  X_initial <- t(model$X_initial)
  X_transition <- aperm(model$X_transition, c(3, 1, 2))
  X_emission <- aperm(model$X_emission, c(3, 1, 2))
  X_cluster <- t(model$X_cluster)
  theta_raw <- model$coefficients$theta_raw
  initial_probs <- vector("list", D)
  transition_probs <- vector("list", D)
  emission_probs <- vector("list", D)
  for (d in seq_len(D)) {
    beta_i_raw <- stan_to_cpp_initial(
      matrix(
        model$coefficients$beta_i_raw[d, ,], 
        S - 1, nrow(X_initial)
      )
    )
    beta_s_raw <- stan_to_cpp_transition(
      array(
        model$coefficients$beta_s_raw[d, , ,], 
        dim = c(S, S - 1, nrow(X_transition))
      )
    )
    beta_o_raw <- stan_to_cpp_emission(
      if (C == 1) {
        array(
          model$coefficients$beta_o_raw[d, , ,],
          dim = c(S, M - 1, nrow(X_emission))
        )
      } else {
        model$coefficients$beta_o_raw[d, ]
      },
      1,
      C > 1
    )
    initial_probs[[d]] <- get_pi(beta_i_raw, X_initial, 0)
    transition_probs[[d]] <- get_A(beta_s_raw, X_transition, 0)
    emission_probs[[d]] <- if (C == 1) {
      get_B(beta_o_raw, X_emission, 0, 0) 
    } else {
      get_multichannel_B(beta_o_raw, X_emission, S, C, M, 0, 0) 
    }
  }
  if (C == 1) {
    ids <- rownames(model$observations)
    times <- colnames(model$observations)
    symbol_names <- list(model$symbol_names)
  } else {
    ids <- rownames(model$observations[[1]])
    times <- colnames(model$observations[[1]])
    symbol_names <- model$symbol_names
  }
  initial_probs <- data.frame(
    cluster = rep(model$cluster_names, each = S * N),
    id = rep(ids, each = S),
    state = model$state_names,
    estimate = unlist(initial_probs)
  )
  colnames(initial_probs)[2] <- model$id_variable
  transition_probs <- data.frame(
    cluster = rep(model$cluster_names, each = S^2 * T * N),
    id = rep(ids, each = S^2 * T),
    time = rep(times, each = S^2),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = unlist(transition_probs)
  )
  colnames(transition_probs)[2] <- model$id_variable
  colnames(transition_probs)[3] <- model$time_variable
  emission_probs <- data.frame(
    cluster = rep(model$cluster_names, each = S * sum(M) * T * N),
    id = unlist(lapply(seq_len(C), function(i) rep(ids, each = S * M[i] * T))),
    time = unlist(lapply(seq_len(C), function(i) rep(times, each = S * M[i]))),
    state = model$state_names,
    channel = rep(model$channel_names, S * M * T * N),
    observation = rep(unlist(symbol_names), each = S),
    estimate = unlist(emission_probs)
  )
  colnames(emission_probs)[2] <- model$id_variable
  colnames(emission_probs)[3] <- model$time_variable
  if (C == 1) emission_probs$channel <- NULL
  cluster_probs <- data.frame(
    cluster = model$cluster_names,
    id = rep(ids, each = D),
    estimate = c(get_omega(theta_raw, X_cluster, 0))
  )
  
  if (nsim > 0) {
    out <- sample_parameters(model, nsim, probs)
    for(i in seq_along(probs)) {
      initial_probs[paste0("q", 100 * probs[i])] <- out$quantiles_pi[, i]
    }
    for(i in seq_along(probs)) {
      transition_probs[paste0("q", 100 * probs[i])] <- out$quantiles_A[, i]
    }
    for(i in seq_along(probs)) {
      emission_probs[paste0("q", 100 * probs[i])] <- out$quantiles_B[, i]
    }
    for(i in seq_along(probs)) {
      cluster_probs[paste0("q", 100 * probs[i])] <- out$quantiles_omega[, i]
    }
  }
  rownames(initial_probs) <- NULL
  rownames(transition_probs) <- NULL
  rownames(emission_probs) <- NULL
  rownames(cluster_probs) <- NULL
  list(
    initial_probs = initial_probs, 
    transition_probs = remove_voids(model, transition_probs),
    emission_probs = remove_voids(model, emission_probs),
    cluster_probs = cluster_probs
  )
}
