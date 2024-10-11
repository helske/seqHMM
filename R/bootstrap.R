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
#' Bootstrap Sampling of NHMM Coefficients
#' 
#' @param model An `nhmm` or `mnhmm` object.
#' @param B number of bootstrap samples.
#' @param method Either `"nonparametric"` or `"parametric"`, to define whether 
#' nonparametric or parametric bootstrap should be used. The former samples 
#' sequences with replacement, whereas the latter simulates new datasets based 
#' on the model.
#' @param penalty penalty term for model estimation. By default, same penalty is used 
#' as was in model estimation by [estimate_nhmm()] or [estimate_mnhmm()].
#' @param verbose Should the progress bar be displayed? Default is `FALSE`.
#' @param ... Additional arguments to [nloptr()].
#' @return The original model with additional element `model$boot`, which 
#' contains A n * B matrix where n is the number of coefficients. The order of 
#' coefficients match `unlist(model$gammas)`. 
#' @rdname bootstrap
#' @export
bootstrap_coefs <- function(model, ...) {
  UseMethod("bootstrap_coefs", model)
}
#' @rdname bootstrap
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
  gamma_pi <- replicate(B, gammas_mle$pi, simplify = FALSE)
  gamma_A <- replicate(B, gammas_mle$A, simplify = FALSE)
  gamma_B <- replicate(B, gammas_mle$B, simplify = FALSE)
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
  if (method == "nonparametric") {
    for (i in seq_len(B)) {
      mod <- bootstrap_model(model)
      fit <- fit_nhmm(mod, init, 0, 0, 1, penalty, ...)
      fit$gammas <- permute_states(fit$gammas, gammas_mle)
      gamma_pi[[i]] <- fit$gammas$pi
      gamma_A[[i]] <- fit$gammas$A
      gamma_B[[i]] <- fit$gammas$B
      if (verbose) utils::setTxtProgressBar(pb, 100 * i/B)
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
      fit$gammas <- permute_states(fit$gammas, gammas_mle)
      gamma_pi[[i]] <- fit$gammas$pi
      gamma_A[[i]] <- fit$gammas$A
      gamma_B[[i]] <- fit$gammas$B
      if (verbose) utils::setTxtProgressBar(pb, 100 * i/B)
    }
  }
  if (verbose) close(pb)
  model$boot <- list(gamma_pi = gamma_pi, gamma_A = gamma_A, gamma_B = gamma_B)
  model
}
#' @rdname bootstrap
#' @export
bootstrap_coefs.mnhmm <- function(model, B = 1000, 
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
  gamma_pi <- replicate(B, gammas_mle$pi, simplify = FALSE)
  gamma_A <- replicate(B, gammas_mle$A, simplify = FALSE)
  gamma_B <- replicate(B, gammas_mle$B, simplify = FALSE)
  gamma_omega <- replicate(B, gammas_mle$omega, simplify = FALSE)
  D <- model$n_clusters
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
  if (method == "nonparametric") {
    for (i in seq_len(B)) {
      mod <- bootstrap_model(model)
      fit <- fit_mnhmm(mod, init, 0, 0, 1, penalty, ...)
      for (j in seq_len(D)) {
        out <- permute_states(
          lapply(fit$gammas[-4], "[[", j), 
          lapply(gammas_mle[-4], "[[", j)
        )
        fit$gammas$pi[[j]] <- out$pi
        fit$gammas$A[[j]] <- out$A
        fit$gammas$B[[j]] <- out$B
      }
      gamma_pi[[i]] <- fit$gammas$pi
      gamma_A[[i]] <- fit$gammas$A
      gamma_B[[i]] <- fit$gammas$B
      gamma_omega[[i]] <- fit$gammas$omega
      if (verbose) utils::setTxtProgressBar(pb, 100 * i/B)
    }
  } else {
    N <- model$n_sequences
    T_ <- model$sequence_lengths
    M <- model$n_symbols
    S <- model$n_states
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
      fit <- fit_mnhmm(mod, init, 0, 0, 1, penalty, ...)
      for (j in seq_len(D)) {
        out <- permute_states(
          lapply(fit$gammas[-4], "[[", j), 
          lapply(gammas_mle[-4], "[[", j)
        )
        fit$gammas$pi[[j]] <- out$pi
        fit$gammas$A[[j]] <- out$A
        fit$gammas$B[[j]] <- out$B
      }
      gamma_pi[[i]] <- fit$gammas$pi
      gamma_A[[i]] <- fit$gammas$A
      gamma_B[[i]] <- fit$gammas$B
      gamma_omega[[i]] <- fit$gammas$omega
      if (verbose) utils::setTxtProgressBar(pb, 100 * i/B)
    }
  }
  if (verbose) close(pb)
  model$boot <- list(
    gamma_pi = gamma_pi, gamma_A = gamma_A, gamma_B = gamma_B,
    gamma_omega = gamma_omega
  )
  model
}
