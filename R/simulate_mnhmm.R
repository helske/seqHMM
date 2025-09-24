#' Simulate Mixture Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' mixture non-homogeneous hidden Markov model.
#'
#' @param n_clusters The number of clusters/mixtures.
#' @param cluster_formula of class [formula()] for the mixture probabilities.
#' @param coefs Same as argument `inits` in [estimate_mnhmm()]. If `NULL`, 
#' (default), the model parameters are generated randomly. If you want to 
#' simulate new sequences based on an estimated model `fit`, you can use 
#' `coefs = fit$etas` and `init_sd = 0`.
#' @param init_sd Standard deviation of the normal distribution used to 
#' generate random coefficients. Default is `2` when `coefs` is `NULL` and `0` 
#' otherwise.
#' @inheritParams simulate_nhmm
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_mnhmm <- function(
    n_states, n_clusters, emission_formula, initial_formula = ~1, 
    transition_formula = ~1, cluster_formula = ~ 1, data, id, time, 
    coefs = NULL, init_sd = 2 * is.null(coefs), check_rank = NULL) {
  force(init_sd)
  stopifnot_(
    !missing(data),
    "{.arg data} is missing, use {.fn simulate_mhmm} instead."
  )
  stopifnot_(
    !missing(emission_formula),
    "Argument {.arg emission_formula} is missing."
  )
  if (inherits(emission_formula, "formula")) {
    responses <- get_responses(emission_formula)
    C <- length(responses)
    if (C > 1L) {
      rhs <- deparse1(emission_formula[[3L]])
      emission_formula <- lapply(
        responses, \(y) stats::as.formula(
          paste(y, " ~ ", rhs), 
          env = environment(emission_formula)
        )
      )
    } else {
      emission_formula <- list(emission_formula)
    }
  } else {
    responses <- vapply(emission_formula, get_responses, allow_mv = FALSE, "")
    C <- length(responses)
  }
  data <- .check_data(data, id, time, responses)
  for (y in responses) {
    l <- as.factor(levels(data[[y]]))
    data[, y := fifelse(is.na(y[1]), l[1], y[1]), by = id, env = list(y = y), 
         showProgress = FALSE]
  }
  cluster_names <- paste("Cluster", seq_len(n_clusters))
  if (is.null(coefs)) {
    coefs <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL, 
      cluster_probs = NULL)
  } else {
    if (is.null(coefs$initial_probs)) coefs$initial_probs <- NULL
    if (is.null(coefs$transition_probs)) coefs$transition_probs <- NULL
    if (is.null(coefs$emission_probs)) coefs$emission_probs <- NULL
    if (is.null(coefs$cluster_probs)) coefs$cluster_probs <- NULL
  }
  model <- build_mnhmm(
    n_states, n_clusters, emission_formula, initial_formula, transition_formula, 
    cluster_formula, data, id, time, coefs = coefs, scale = FALSE, 
    check = check_rank, drop_levels = FALSE
  )
  model$etas <- create_initial_values(coefs, model, init_sd)
  model$gammas$gamma_pi <- drop(eta_to_gamma_mat_field(model$etas$eta_pi))
  model$gammas$gamma_A <- drop(eta_to_gamma_cube_field(model$etas$eta_A))
  model$gammas$gamma_B <- split(
    eta_to_gamma_cube_2d_field(model$etas$eta_B), 
    seq_len(n_clusters)
  )
  model$gammas$gamma_omega <- eta_to_gamma_mat(model$etas$eta_omega)
  
  if (inherits(model, "fanhmm")) {
    W <- update_W_for_fanhmm(model)
    out <- Rcpp_simulate_mfanhmm(
      create_obs(model), model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$gammas$gamma_omega, model$prior_obs, model$W_X_B, W$W_A, W$W_B
    )
    if (identical(model$prior_obs, 0L)) {
      .idx <- setdiff(
        seq_row(data), 
        cumsum(c(1, head(model$sequence_lengths + 1L, -1)))
      ) 
    } else {
      .idx <- NULL
    }
    for (i in seq_len(model$n_channels)) {
      y <- unlist(lapply(out$observations, \(y) y[i, ] + 1L))
      set(data, i = .idx, j = model$responses[i], 
          value = model$symbol_names[[i]][y])
    }
  } else {
    out <- Rcpp_simulate_mnhmm(
      create_obs(model), model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, model$gammas$gamma_omega
    )
    for (i in seq_len(model$n_channels)) {
      y <- unlist(lapply(out$observations, \(y) y[i, ] + 1L))
      set(data, j = model$responses[i], value = model$symbol_names[[i]][y])
    }
  }
  state_names <- paste0(
    rep(model$cluster_names, each = model$n_states), ": ",
    unlist(model$state_names)
  )
  states <- cbind(
    model$data[, list(id, time), env = list(id = id, time = time)], 
    state = state_names[unlist(out$states) + 1L]
  )
  attr(model$X_omega, "X_mean") <- TRUE
  attr(model$X_omega, "R_inv") <- NULL
  attr(model$X_pi, "X_mean") <- TRUE
  attr(model$X_pi, "R_inv") <- NULL
  attr(model$X_A, "X_mean") <- TRUE
  attr(model$X_A, "R_inv") <- NULL
  for (i in seq_along(responses)) {
    attr(model$X_B[[i]], "X_mean") <- TRUE
    attr(model$X_B[[i]], "R_inv") <- NULL
  }
  model <- update(model, data)
  tQs <- t(create_Q(n_states))
  if (!attr(model$X_pi, "icpt_only")) {
    model$gammas$gamma_pi <- gamma_to_gamma_std(
      model$gammas$gamma_pi, solve(attr(model$X_pi, "R_inv")), 
      attr(model$X_pi, "coef_names"), attr(model$X_pi, "X_mean")
    )
    for (d in seq_len(n_clusters)) {
      model$etas$eta_pi[[d]][] <- tQd %*% model$gammas$gamma_pi[[d]]
    }
  }
  if (!attr(model$X_A, "icpt_only")) {
    model$gammas$gamma_A <- gamma_to_gamma_std(
      model$gammas$gamma_A, solve(attr(model$X_A, "R_inv")), 
      attr(model$X_A, "coef_names"), attr(model$X_A, "X_mean")
    )
    for (d in seq_len(n_clusters)) {
      for (s in seq_len(n_states)) {
        model$etas$eta_A[[d]][, , s] <- tQs %*% model$gammas$gamma_A[[d]][, , s]
      }
    }
  }
  for (i in seq_along(responses)) {
    if (!attr(model$X_B[[i]], "icpt_only")) {
      tQm <- t(create_Q(model$n_symbols[i]))
      coef_names <- attr(model$X_B[[i]], "coef_names")
      R_inv <- solve(attr(model$X_B[[i]], "R_inv"))
      X_mean <- attr(model$X_B[[i]], "X_mean")
      coef_names <- attr(model$X_B[[i]], "coef_names")
      for (d in seq_len(n_clusters)) {
        model$gammas$gamma_B[[d]][[i]] <- gamma_to_gamma_std(
          model$gammas$gamma_B[[d]][[i]], R_inv, coef_names, X_mean
        )
        for (s in seq_len(n_states)) {
          model$etas$eta_B[[d]][[i]][, , s] <- tQm %*% model$gammas$gamma_B[[d]][[i]][, , s]
        }
      }
    }
  }
  if (!attr(model$X_omega, "icpt_only")) {
    tQd <- t(create_Q(n_clusters))
    model$gammas$gamma_omega <- gamma_to_gamma_std(
      model$gammas$gamma_omega, solve(attr(model$X_omega, "R_inv")), 
      attr(model$X_omega, "coef_names"), attr(model$X_omega, "X_mean")
    )
    model$etas$eta_omega[] <- tQd %*% model$gammas$gamma_omega
  }
  list(model = model, states = states, data = data)
}
