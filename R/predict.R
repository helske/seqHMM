#' Predict method for non-homogeneous hidden Markov models
#' 
#' @param object A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param newdata Optional data frame which is used for prediction.
#' @param nsim Non-negative integer defining the number of samples from the
#' normal approximation of the model coefficients.
#' @param probs Vector defining the quantiles of interest.
#' @param return_samples Logical indicating whether to return samples or
#' quantiles. Default is `FALSE`.
#' @param ... Ignored.
#' @export
predict.nhmm <- function(
    object, newdata = NULL, nsim = 0, 
    probs = c(0.025, 0.5, 0.975), return_samples = FALSE, ...) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  if (!is.null(newdata)) {
    time <- object$time_variable
    id <- object$id_variable
    stopifnot_(
      is.data.frame(newdata), 
      "Argument {.arg newdata} must be a {.cls data.frame} object."
    )
    stopifnot_(
      !is.null(newdata[[id]]), 
      "Can't find grouping variable {.var {id}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[time]]), 
      "Can't find time index variable {.var {time}} in {.arg newdata}."
    )
    object <- update(object, newdata = newdata)
  }
  S <- object$n_states
  M <- object$n_symbols
  C <- object$n_channels
  N <- object$n_sequences
  T_ <- object$length_of_sequences
  beta_i_raw <- stan_to_cpp_initial(
    object$coefficients$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    object$coefficients$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    object$coefficients$beta_o_raw,
    1,
    C > 1
  )
  X_initial <- t(object$X_initial)
  X_transition <- aperm(object$X_transition, c(3, 1, 2))
  X_emission <- aperm(object$X_emission, c(3, 1, 2))
  initial_probs <- get_pi(beta_i_raw, X_initial, 0)
  transition_probs <- get_A(beta_s_raw, X_transition, 0)
  emission_probs <- if (C == 1) {
    get_B(beta_o_raw, X_emission, 0, 0) 
  } else {
    get_multichannel_B(beta_o_raw, X_emission, S, C, M, 0, 0) 
  } 
  if (C == 1) {
    ids <- rownames(object$observations)
    times <- colnames(object$observations)
    symbol_names <- list(object$symbol_names)
  } else {
    ids <- rownames(object$observations[[1]])
    times <- colnames(object$observations[[1]])
    symbol_names <- object$symbol_names
  }
  out <- list()
  out$initial_probs <- data.frame(
    id = rep(ids, each = S),
    state = object$state_names,
    estimate = c(initial_probs)
  )
  colnames(out$initial_probs)[1] <- object$id_variable
  out$transition_probs <- data.frame(
    id = rep(ids, each = S^2 * T_),
    time = rep(times, each = S^2),
    state_from = object$state_names,
    state_to = rep(object$state_names, each = S),
    estimate = unlist(transition_probs)
  )
  colnames(out$transition_probs)[1] <- object$id_variable
  colnames(out$transition_probs)[2] <- object$time_variable
  out$emission_probs <- do.call(
    "rbind", 
    lapply(seq_len(C), function(i) {
      data.frame(
        id = rep(ids, each = S * M[i] * T_),
        time = rep(times, each = S * M[i]),
        state = object$state_names,
        channel = object$channel_names[i],
        observation = rep(symbol_names[[i]], each = S),
        estimate = unlist(emission_probs[((i - 1) * N + 1):(i * N)])
      )
    })
  )
  colnames(out$emission_probs)[1] <- object$id_variable
  colnames(out$emission_probs)[2] <- object$time_variable 
  if (C == 1) out$emission_probs$channel <- NULL
  
  if (nsim > 0) {
    samples <- sample_parameters_nhmm(object, nsim)
    if (return_samples) {
      out$samples <- samples_to_df(object, samples)
    } else {
      out$quantiles <-  list(
        initial_probs = fast_quantiles(samples$pi, probs),
        transition_probs = fast_quantiles(samples$A, probs),
        emission_probs = fast_quantiles(samples$B, probs)
      )
    }
  }
  out
}
#' @export
predict.mnhmm <- function(
    object, newdata = NULL, nsim = 0, 
    probs = c(0.025, 0.5, 0.975), return_samples = FALSE, ...) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  if (!is.null(newdata)) {
    time <- object$time_variable
    id <- object$id_variable
    stopifnot_(
      is.data.frame(newdata), 
      "Argument {.arg newdata} must be a {.cls data.frame} object."
    )
    stopifnot_(
      !is.null(newdata[[id]]), 
      "Can't find grouping variable {.var {id}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[time]]), 
      "Can't find time index variable {.var {time}} in {.arg newdata}."
    )
    object <- update(object, newdata = newdata)
  }
  T_ <- object$length_of_sequences
  N <- object$n_sequences
  S <- object$n_states
  M <- object$n_symbols
  C <- object$n_channels
  D <- object$n_clusters
  X_initial <- t(object$X_initial)
  X_transition <- aperm(object$X_transition, c(3, 1, 2))
  X_emission <- aperm(object$X_emission, c(3, 1, 2))
  X_cluster <- t(object$X_cluster)
  theta_raw <- object$coefficients$theta_raw
  initial_probs <- vector("list", D)
  transition_probs <- vector("list", D)
  emission_probs <- vector("list", D)
  for (d in seq_len(D)) {
    beta_i_raw <- stan_to_cpp_initial(
      matrix(
        object$coefficients$beta_i_raw[d, ,], 
        S - 1, nrow(X_initial)
      )
    )
    beta_s_raw <- stan_to_cpp_transition(
      array(
        object$coefficients$beta_s_raw[d, , ,], 
        dim = c(S, S - 1, nrow(X_transition))
      )
    )
    beta_o_raw <- stan_to_cpp_emission(
      if (C == 1) {
        array(
          object$coefficients$beta_o_raw[d, , ,],
          dim = c(S, M - 1, nrow(X_emission))
        )
      } else {
        object$coefficients$beta_o_raw[d, ]
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
    ids <- rownames(object$observations)
    times <- colnames(object$observations)
    symbol_names <- list(object$symbol_names)
  } else {
    ids <- rownames(object$observations[[1]])
    times <- colnames(object$observations[[1]])
    symbol_names <- object$symbol_names
  }
  out <- list()
  out$initial_probs <- data.frame(
    cluster = rep(object$cluster_names, each = S * N),
    id = rep(ids, each = S),
    state = object$state_names,
    estimate = unlist(initial_probs)
  )
  colnames(out$initial_probs)[2] <- object$id_variable
  out$transition_probs <- data.frame(
    cluster = rep(object$cluster_names, each = S^2 * T_ * N),
    id = rep(ids, each = S^2 * T_),
    time = rep(times, each = S^2),
    state_from = object$state_names,
    state_to = rep(object$state_names, each = S),
    estimate = unlist(transition_probs)
  )
  colnames(out$transition_probs)[2] <- object$id_variable
  colnames(out$transition_probs)[3] <- object$time_variable
  out$emission_probs <- data.frame(
    cluster = rep(object$cluster_names, each = S * sum(M) * T_ * N),
    id = unlist(lapply(seq_len(C), function(i) rep(ids, each = S * M[i] * T_))),
    time = unlist(lapply(seq_len(C), function(i) rep(times, each = S * M[i]))),
    state = object$state_names,
    channel = rep(object$channel_names, S * M * T_ * N),
    observation = rep(unlist(symbol_names), each = S),
    estimate = unlist(emission_probs)
  )
  colnames(out$emission_probs)[2] <- object$id_variable
  colnames(out$emission_probs)[3] <- object$time_variable
  if (C == 1) emission_probs$channel <- NULL
  out$cluster_probs <- data.frame(
    cluster = object$cluster_names,
    id = rep(ids, each = D),
    estimate = c(get_omega(theta_raw, X_cluster, 0))
  )
  if (nsim > 0) {
    samples <- sample_parameters_mnhmm(object, nsim)
    if (return_samples) {
      out$samples <- samples_to_df(object, samples)
    } else {
      out$quantiles <-  list(
        initial_probs = fast_quantiles(samples$pi, probs),
        transition_probs = fast_quantiles(samples$A, probs),
        emission_probs = fast_quantiles(samples$B, probs),
        cluster_probs = fast_quantiles(samples$omega, probs)
      )
    }
  }
  out
}
