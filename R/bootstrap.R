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
  model$sequence_lengths <- model$sequence_lengths[idx]
  model
}
permute_states <- function(gammas_boot, gammas_mle) {
  C <- if(is.list(gammas_mle$B)) length(gammas_mle$B) else 1
  if (C == 1) {
    m <- cost_matrix_singlechannel(
      gammas_boot$pi, gammas_mle$pi,
      gammas_boot$A, gammas_mle$A,
      gammas_boot$B, gammas_mle$B
    )
  } else {
    m <- cost_matrix_multichannel(
      gammas_boot$pi, gammas_mle$pi,
      gammas_boot$A, gammas_mle$A,
      gammas_boot$B, gammas_mle$B
    )
  }
  perm <- RcppHungarian::HungarianSolver(m)$pairs[, 2]
  gammas_boot$pi <- gammas_boot$pi[perm, , drop = FALSE]
  gammas_boot$A <- gammas_boot$A[perm, , perm, drop = FALSE]
  if (C == 1) {
    gammas_boot$B <- gammas_boot$B[, , perm, drop = FALSE]
  } else {
    for (c in seq_len(C)) {
      gammas_boot$B[[c]] <- gammas_boot$B[[c]][, , perm, drop = FALSE]
    }
  }
  gammas_boot
}
#' @export
bootstrap_coefs.nhmm <- function(model, B = 1000, 
                                 method = c("nonparametric", "parametric"),
                                 penalty, verbose = FALSE, ...) {
  method <- match.arg(method)
  stopifnot_(
    checkmate::test_int(x = B, lower = 0L), 
    "Argument {.arg B} must be a single positive integer."
  )
  init <- model$etas
  if (missing(penalty)) {
    penalty <- model$estimation_results$penalty
  }
  gammas_mle <- model$gammas
  
  coefs <- matrix(NA, length(unlist(gammas_mle)), B)
  if (method == "nonparametric") {
    for (i in seq_len(B)) {
      mod <- bootstrap_model(model)
      fit <- fit_nhmm(mod, init, 0, 0, 1, penalty, ...)
      coefs[, i] <- unlist(permute_states(fit$gammas, gammas_mle))
      if(verbose) print(paste0("Bootstrap replication ", i, " complete."))
    }
  } else {
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
      fit <- fit_nhmm(mod, init, 0, 0, 1, penalty, ...)
      coefs[, i] <- unlist(permute_states(fit$gammas, gammas_mle))
      print(paste0("Bootstrap replication ", i, " complete."))
    }
  }
  return(coefs)
}
#' @export
bootstrap_coefs.mnhmm <- function(model, B = 1000, 
                                  method = c("nonparametric", "parametric"),
                                  verbose = FALSE) {
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
      print(paste0("Bootstrap replication ", i, " complete."))
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
      print(paste0("Bootstrap replication ", i, " complete."))
    }
  }
  return(coefs)
}
