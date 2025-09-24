#' Simulate Non-homogeneous Hidden Markov Models
#'
#' Simulate sequences of observed and hidden states given the parameters of a 
#' non-homogeneous hidden Markov model.
#' 
#' @param data A data frame containing the variables used in the model 
#' formulas. Note that this should also include also the response variable(s), 
#' which are used to define the number of observed symbols (using [levels()]) 
#' and the length of sequences. The actual values of the response variables 
#' does not matter though, as they are replaced by the simulated values. The 
#' exception is the first time point in FAN-HMM case: If the `emission_formula` 
#' contains lagged responses, the response variable values at the first time 
#' point are used to define the emissions at the second time point, and the 
#' simulations are done from the second time point onward. This matches the 
#' case `prior_obs = "fixed"` in [estimate_nhmm()]. Note that compared to 
#' `estimate_*` functions, unused factor levels are not automatically dropped 
#' from `data`.
#' @param coefs Same as argument `inits` in [estimate_nhmm()]. If `NULL`,
#' (default), the model parameters are generated randomly. If you want to 
#' simulate new sequences based on an estimated model `fit`, you can use 
#' `coefs = fit$etas` and `init_sd = 0`.
#' @param init_sd Standard deviation of the normal distribution used to 
#' generate random coefficients. Default is `2` when `coefs` is `NULL` and `0` 
#' otherwise.
#' @inheritParams estimate_nhmm
#'
#' @return A list with the model used in simulation as well as the simulated 
#' hidden state sequences.
#' @export
simulate_nhmm <- function(
    n_states, emission_formula, initial_formula = ~1, transition_formula = ~1, 
    data, id, time, coefs = NULL, init_sd = 2 * is.null(coefs), 
    check_rank = NULL) {
  
  force(init_sd)
  stopifnot_(
    !missing(data),
    "{.arg data} is missing, use {.fn simulate_hmm} instead."
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
    l <- as.factor(levels(data[[y]]))[1]
    data[, y := fifelse(is.na(y[1]), l, y[1]), by = id, env = list(y = y), 
         showProgress = FALSE]
  }
  if (is.null(coefs)) {
    coefs <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL)
  } else {
    if (is.null(coefs$initial_probs)) coefs$initial_probs <- NULL
    if (is.null(coefs$transition_probs)) coefs$transition_probs <- NULL
    if (is.null(coefs$emission_probs)) coefs$emission_probs <- NULL
  }
  model <- build_nhmm(
    n_states, emission_formula, initial_formula, transition_formula, 
    data, id, time, coefs = coefs, scale = FALSE, 
    check = check_rank, drop_levels = FALSE
  )
  model$etas <- create_initial_values(coefs, model, init_sd)
  model$gammas$gamma_pi <- eta_to_gamma_mat(model$etas$eta_pi)
  model$gammas$gamma_A <- eta_to_gamma_cube(model$etas$eta_A)
  model$gammas$gamma_B <- drop(eta_to_gamma_cube_field(model$etas$eta_B))
  if (inherits(model, "fanhmm")) {
    W <- update_W_for_fanhmm(model)
    out <- Rcpp_simulate_fanhmm(
      create_obs(model), model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$prior_obs, model$W_X_B, W$W_A, W$W_B
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
    out <- Rcpp_simulate_nhmm(
      create_obs(model), model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B
    )
    for (i in seq_len(model$n_channels)) {
      y <- unlist(lapply(out$observations, \(y) y[i, ] + 1L))
      set(data, j = model$responses[i], value = model$symbol_names[[i]][y])
    }
  }
  states <- cbind(
    model$data[, list(id, time), env = list(id = id, time = time)], 
    state = model$state_names[unlist(out$states) + 1L]
  )
  
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
    model$etas$eta_pi[] <- tQs %*% model$gammas$gamma_pi
  }
  if (!attr(model$X_A, "icpt_only")) {
    model$gammas$gamma_A <- gamma_to_gamma_std(
      model$gammas$gamma_A, solve(attr(model$X_A, "R_inv")), 
      attr(model$X_A, "coef_names"), attr(model$X_A, "X_mean")
    )
    for (s in seq_len(n_states)) {
      model$etas$eta_A[, , s] <- tQs %*% model$gammas$gamma_A[, , s]
    }
  }
  for (i in seq_along(responses)) {
    if (!attr(model$X_B[[i]], "icpt_only")) {
      tQm <- t(create_Q(model$n_symbols[i]))
      model$gammas$gamma_B[[i]] <- gamma_to_gamma_std(
        model$gammas$gamma_B[[i]], solve(attr(model$X_B[[i]], "R_inv")), 
        attr(model$X_B[[i]], "coef_names"), attr(model$X_B[[i]], "X_mean")
      )
      for (s in seq_len(n_states)) {
        model$etas$eta_B[[i]][, , s] <- tQm %*% model$gammas$gamma_B[[i]][, , s]
      }
    }
  }
  list(model = model, states = states, data = data)
}
