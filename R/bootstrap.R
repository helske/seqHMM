#' Bootstrap observations
#' @noRd
bootstrap_model <- function(model) {
  idx <- sample.int(model$n_sequences, replace = TRUE)
  if (model$n_channels  == 1) {
    model$observations[] <- model$observations[idx, ]
  } else {
    for(i in seq_len(model$n_channels)) {
      model$observations[[i]][] <- model$observations[[i]][idx, ,]
    }
  }
  model$X_pi[] <- model$X_pi[, idx]
  model$X_A[] <- model$X_A[, , idx]
  model$X_B[] <- model$X_B[, , idx]
  if (!is.null(model$X_omega)) {
    model$X_omega[] <- model$X_omega[, idx]
  }
  model$sequence_lengths <- model$sequence_lengths[idx]
  model
}
#' Permute states of bootstrap sample to match MLE
#' @noRd
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
#' Permute clusters of bootstrap sample to match MLE
#' @noRd
permute_clusters <- function(model, pcp_mle) {
  pcp <- posterior_cluster_probabilities(model)
  m <- cost_matrix_clusters(pcp, pcp_mle)
  perm <- RcppHungarian::HungarianSolver(m)$pairs[, 2]
  model$gammas$omega[perm, , drop = FALSE]
  model$gammas$pi <- model$gammas$pi[perm]
  model$gammas$A <- model$gammas$A[perm]
  model$gammas$B <- model$gammas$B[perm]
  model
}
#' Bootstrap Sampling of NHMM Coefficients
#' 
#' It is possible to parallelize the bootstrap runs using the `future` package, 
#' e.g., by calling `future::plan(multisession, workers = 2)` before 
#' `bootstrap_coefs()`. See [future::plan()] for details.
#' 
#' `bootstrap_coefs()` is compatible with `progressr` package, so you can use 
#' `progressr::with_progress(bootstrap_coefs(fit))` to track the progress of 
#' bootstrapping.
#' 
#' @param model An `nhmm` or `mnhmm` object.
#' @param B number of bootstrap samples.
#' @param type Either `"nonparametric"` or `"parametric"`, to define whether 
#' nonparametric or parametric bootstrap should be used. The former samples 
#' sequences with replacement, whereas the latter simulates new datasets based 
#' on the model.
#' @param method Estimation method used in bootstrapping. Defaults to 
#' `"EM-LBFGS"`.
#' @param ... Additional arguments to [estimate_nhmm()] or [estimate_mnhmm()].
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
                                 type = c("nonparametric", "parametric"),
                                 method = "EM-LBFGS", ...) {
  type <- match.arg(type)
  stopifnot_(
    checkmate::test_int(x = B, lower = 0L), 
    "Argument {.arg B} must be a single positive integer."
  )
  init <- model$etas
  gammas_mle <- model$gammas
  lambda <- model$estimation_results$lambda
  pseudocount <- model$estimation_results$pseudocount
  p <- progressr::progressor(along = seq_len(B))
  if (type == "nonparametric") {
    out <- future.apply::future_lapply(
      seq_len(B), function(i) {
        mod <- bootstrap_model(model)
        fit <- fit_nhmm(mod, init, init_sd = 0, restarts = 0, lambda = lambda, 
                        method = method, pseudocount = pseudocount, ...)
        fit$gammas <- permute_states(fit$gammas, gammas_mle)
        p()
        fit$gammas
      }, future.seed = TRUE
    )
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
    out <- future.apply::future_lapply(
      seq_len(B), function(i) {
        mod <- simulate_nhmm(
          N, T_, M, S, formula_pi, formula_A, formula_B,
          data = d, time, id, init)$model
        fit <- fit_nhmm(mod, init, init_sd = 0, restarts = 0, lambda = lambda, 
                        method = method, pseudocount = pseudocount, ...)
        fit$gammas <- permute_states(fit$gammas, gammas_mle)
        p()
        fit$gammas
      }, future.seed = TRUE
    )
  }
  model$boot <- list(
    gamma_pi = lapply(out, "[[", "pi"), 
    gamma_A = lapply(out, "[[", "A"), 
    gamma_B = lapply(out, "[[", "B")
  )
  model
}
#' @rdname bootstrap
#' @export
bootstrap_coefs.mnhmm <- function(model, B = 1000, 
                                  type = c("nonparametric", "parametric"),
                                  method = "EM-LBFGS", ...) {
  type <- match.arg(type)
  stopifnot_(
    checkmate::test_int(x = B, lower = 0L), 
    "Argument {.arg B} must be a single positive integer."
  )
  init <- model$etas
  gammas_mle <- model$gammas
  pcp_mle <- posterior_cluster_probabilities(model)
  lambda <- model$estimation_results$lambda
  pseudocount <- model$estimation_results$pseudocount
  D <- model$n_clusters
  p <- progressr::progressor(along = seq_len(B))
  if (type == "nonparametric") {
    out <- future.apply::future_lapply(
      seq_len(B), function(i) {
        mod <- bootstrap_model(model)
        fit <- fit_mnhmm(mod, init, init_sd = 0, restarts = 0, lambda = lambda, 
                         method = method, pseudocount = pseudocount, ...)
        fit <- permute_clusters(fit, pcp_mle)
        for (j in seq_len(D)) {
          out <- permute_states(
            lapply(fit$gammas[c("pi", "A", "B")], "[[", j), 
            lapply(gammas_mle[c("pi", "A", "B")], "[[", j)
          )
          fit$gammas$pi[[j]] <- out$pi
          fit$gammas$A[[j]] <- out$A
          fit$gammas$B[[j]] <- out$B
        }
        p()
        fit$gammas
      }, future.seed = TRUE
    )
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
    out <- future.apply::future_lapply(
      seq_len(B), function(i) {
        mod <- simulate_mnhmm(
          N, T_, M, S, D, formula_pi, formula_A, formula_B, formula_omega,
          data = d, time, id, init)$model
        fit <- fit_mnhmm(mod, init, init_sd = 0, restarts = 0, lambda = lambda, 
                         method = method, pseudocount = pseudocount, ...)
        fit <- permute_clusters(fit, pcp_mle)
        for (j in seq_len(D)) {
          out <- permute_states(
            lapply(fit$gammas[c("pi", "A", "B")], "[[", j), 
            lapply(gammas_mle[c("pi", "A", "B")], "[[", j)
          )
          fit$gammas$pi[[j]] <- out$pi
          fit$gammas$A[[j]] <- out$A
          fit$gammas$B[[j]] <- out$B
        }
        p()
        fit$gammas
      }, future.seed = TRUE
    )
  }
  model$boot <- list(
    gamma_pi = lapply(out, "[[", "pi"), 
    gamma_A = lapply(out, "[[", "A"), 
    gamma_B = lapply(out, "[[", "B"),
    gamma_omega = lapply(out, "[[", "omega")
  )
  model
}
