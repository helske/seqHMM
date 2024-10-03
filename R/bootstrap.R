bootstrap_model <- function(model) {
  idx <- sample.int(model$n_sequences, replace = TRUE)
  if (model$n_channels  == 1) {
    model$observations <- model$observations[idx, , drop = FALSE]
  } else {
    for(i in seq_len(model$n_channels)) {
      model$observations[[i]] <- model$observations[[i]][idx, , drop = FALSE]
    }
  }
  model$X_initial <- model$X_initial[, idx, drop = FALSE]
  model$X_transition <- model$X_transition[, , idx, drop = FALSE]
  model$X_emission <- model$X_emission[, , idx, drop = FALSE]
  if (!is.null(model$X_cluster)) {
    model$X_cluster <- model$X_cluster[, idx, drop = FALSE]
  }
  model
}
#' @export
bootstrap_coefs.nhmm <- function(model, B = 1000, 
                                 method = c("nonparametric", "parametric"),
                                 penalty) {
  method <- match.arg(method)
  stopifnot_(
    checkmate::test_int(x = B, lower = 0L), 
    "Argument {.arg B} must be a single positive integer."
  )
  init <- model$coefficients
  if (missing(penalty)) {
    penalty <- model$estimation_results$penalty
  }
  mle_coefficients <- model$coefficients
  if (method == "nonparametric") {
    coefs <- matrix(NA, length(unlist(init)), B)
    for (i in seq_len(B)) {
      mod <- bootstrap_model(model)
      fit <- fit_nhmm(mod, init, 0, 0, 1, penalty, FALSE)
      m <- cost_matrix(fit$coefficients$eta_pi, mle_coefficients$eta_pi,
                       fit$coefficients$eta_A, mle_coefficients$eta_A,
                       fit$coefficients$eta_B, mle_coefficients$eta_B,
                       fit$X_initial, fit$X_transition, fit$X_emission)
      perm <- RcppHungarian:HungarianSolver(m)$pairs[, 2]
      fit$coefficients$gamma_pi[perm]
      coefs[, i] <- unlist(fit$coefficients)
    }
  } else {
    coefs <- matrix(NA, length(unlist(init)), B)
    N <- model$n_sequences
    T_ <- model$sequence_lengths
    M <- model$n_symbols
    S <- model$n_states
    formula_pi <- model$initial_formula
    formula_A <- model$transition_formula
    formula_B <- model$emission_formula
    d <- model$data
    time <- model$time_variable
    id <- model$id_variable
    for (i in seq_len(B)) {
      mod <- simulate_nhmm(
        N, T_, M, S, formula_pi, formula_A, formula_B,
        data = d, time, id, init)$model
      fit <- fit_nhmm(mod, init, 0, 0, 1, penalty, FALSE)
      coefs[, i] <- unlist(fit$coefficients)
    }
  }
  return(coefs)
}
#' @export
bootstrap_coefs.mnhmm <- function(model, B = 1000, 
                                  method = c("nonparametric", "parametric")) {
  method <- match.arg(method)
  stopifnot_(
    checkmate::test_int(x = B, lower = 0L), 
    "Argument {.arg B} must be a single positive integer."
  )
  init <- model$coefficients
  if (missing(penalty)) {
    penalty <- model$estimation_results$penalty
  }
  if (method == "nonparametric") {
    coefs <- matrix(NA, length(unlist(init)), B)
    for (i in seq_len(B)) {
      mod <- bootstrap_model(model)
      fit <- fit_mnhmm(mod, init, 0, 0, 1, penalty, FALSE)
      coefs[, i] <- unlist(fit$coefficients)
    }
  } else {
    coefs <- matrix(NA, length(unlist(init)), B)
    N <- model$n_sequences
    T_ <- model$sequence_lengths
    M <- model$n_symbols
    S <- model$n_states
    D <- model$n_clusters
    formula_pi <- model$initial_formula
    formula_A <- model$transition_formula
    formula_B <- model$emission_formula
    formula_omega <- model$cluster_formula
    d <- model$data
    time <- model$time_variable
    id <- model$id_variable
    for (i in seq_len(B)) {
      mod <- simulate_mnhmm(
        N, T_, M, S, D, formula_pi, formula_A, formula_B, formula_omega,
        data = d, time, id, init)$model
      fit <- fit_mnhmm(mod, init, 0, 0, 1, penalty, FALSE)
      coefs[, i] <- unlist(fit$coefficients)
    }
  }
  return(coefs)
}
