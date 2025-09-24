#' Bootstrap observations
#' @noRd
bootstrap_model <- function(model, ids) {
  id <- y <- NULL
  idx <- sort(sample.int(model$n_sequences, replace = TRUE))
  model$sequence_lengths <- model$sequence_lengths[idx]
  ids <- ids[idx]
  model$data <- model$data[list(ids), on = model$id_variable, 
                           allow.cartesian = TRUE]
  new_ids <- rep(seq_len(model$n_sequences), times = model$sequence_lengths)
  set(model$data, j = model$id_variable, value = new_ids)
  model$X_pi[] <- model$X_pi[, idx]
  model$X_A[] <- model$X_A[idx]
  for (y in model$responses) {
    model$X_B[[y]][] <- model$X_B[[y]][idx]
  }
  if (inherits(model, "mnhmm")) {
    model$X_omega[] <- model$X_omega[, idx]
  }
  if (inherits(model, "fanhmm") & !identical(model$prior_obs, 0)) {
    for (i in seq_along(model$W_X_B)) {
      for (j in seq_along(model$W_X_B[[1]])) {
        model$W_X_B[[i]][[j]][] <- model$W_X_B[[i]][[j]][idx]
      }
    }
  }
  attr(model, "nobs") <- sum(!is.na(
    model$data[, y, env = list(y = I(model$responses))]
  )) / model$n_channels
  list(model = model, idx = idx)
}
#' Permute the states of NHMM using Hungarian algorithm
#' 
#' This function finds the permutation (according to Hungarian algorithm) 
#' of estimated model coefficients \eqn{\gamma} which best matches some
#' reference values. The main purpose is to permute the bootstrap samples of 
#' model coefficients to match as closely as possible to the maximum likelihood 
#' estimates obtained from the original data.
#' 
#' The cost matrix in Hungarian algorithm is based on the sum of L2 norms of 
#' the differences of estimated and reference values of \eqn{\pi}, \eqn{A} and 
#' \eqn{B}.
#' 
#' This function is mostly meant for internal usage within the bootstrap 
#' methods of `seqHMM`.
#' 
#' @param estimates A list \eqn{\gamma} coefficients as in `model$gammas`,
#' where `model` is an `nhmm` object.
#' @param reference Another list of \eqn{\gamma} coefficients for which to match 
#' the `estimates`.
#' @return Permuted version of `estimates`, with added attribute `permutation` 
#' which contains the permutations used to obtain for example the new 
#' `gamma_pi` as `estimates$gamma_pi[perm, , drop = FALSE]`.
#' @export
permute_states <- function(estimates, reference) {
  C <- length(reference$gamma_B)
  m <- cost_matrix(
    estimates$gamma_pi, reference$gamma_pi, 
    estimates$gamma_A, reference$gamma_A,
    estimates$gamma_B, reference$gamma_B
  )
  perm <- RcppHungarian::HungarianSolver(m)$pairs[, 2]
  estimates$gamma_pi <- estimates$gamma_pi[perm, , drop = FALSE]
  estimates$gamma_A <- estimates$gamma_A[perm, , perm, drop = FALSE]
  for (c in seq_len(C)) {
    estimates$gamma_B[[c]] <- estimates$gamma_B[[c]][, , perm, drop = FALSE]
  }
  attr(estimates, "permutation") <- perm
  estimates
}
#' Permute clusters of bootstrap sample to match MLE
#' This is based on matching the posterior cluster probabilities of the 
#' original sample and model and the ones obtained by using the coefficients 
#' from the bootstrap replicate
#' @noRd
permute_clusters <- function(fit, model, pcp_reference) {
  model$gammas <- fit$gammas
  pcp <- matrix(
    posterior_cluster_probabilities(model)$probability, 
    ncol = model$n_clusters, byrow = TRUE
  )
  m <- cost_matrix_clusters(pcp, pcp_reference)
  perm <- RcppHungarian::HungarianSolver(m)$pairs[, 2]
  fit$gammas$gamma_omega[perm, , drop = FALSE]
  fit$gammas$gamma_pi <- fit$gammas$gamma_pi[perm]
  fit$gammas$gamma_A <- fit$gammas$gamma_A[perm]
  fit$gammas$gamma_B <- fit$gammas$gamma_B[perm]
  fit
}
#' Bootstrap Sampling of NHMM Coefficients
#' 
#' The model estimation for each bootstrap sample uses the same method and 
#' tolerances as the original fit. If you want to change these, you can modify 
#' the elements of the input model such as `model$estimation_results$method` 
#' and `model$controls` before passing it to `bootstrap_coefs()`.
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
#' @param nsim number of bootstrap samples.
#' @param type Either `"nonparametric"` (default) or `"parametric"`, to define 
#' whether nonparametric or parametric bootstrap should be used. The former samples 
#' sequences with replacement, whereas the latter simulates new datasets based 
#' on the model.
#' @param append If `TRUE`, in case the model already contains 
#' bootstrap samples, new samples are appended to `model$boot`. If `FALSE` 
#' (default), old samples are discarded.
#' @param ... Ignored.
#' @return The original model with additional element `model$boot`.
#' @rdname bootstrap
#' @export
bootstrap_coefs <- function(model, ...) {
  UseMethod("bootstrap_coefs", model)
}
#' @rdname bootstrap
#' @export
bootstrap_coefs.nhmm <- function(model, nsim, 
                                 type = c("nonparametric", "parametric"),
                                 append = FALSE, ...) {
  type <- try(match.arg(type), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val nonparametric} or 
    {.val parametric}."
  )
  stopifnot_(
    !missing(nsim) && checkmate::test_int(x = nsim, lower = 0L), 
    "Argument {.arg nsim} must be a single positive integer."
  )
  stopifnot_(
    checkmate::test_logical(x = append), 
    "Argument {.arg append} must be a single logical value."
  )
  init <- model$etas
  gammas_mle <- model$gammas
  lambda <- model$estimation_results$lambda
  bound <- model$estimation_results$bound
  method <- model$estimation_results$method
  p <- progressr::progressor(along = seq_len(nsim))
  original_options <- options(future.globals.maxSize = Inf)
  on.exit(options(original_options))
  control <- model$controls$control
  control$print_level <- 0
  control_mstep <- model$controls$mstep
  control_mstep$print_level <- 0
  if (type == "nonparametric") {
    ids <- unique(model$data[[model$id_variable]])
    out <- future.apply::future_lapply(
      seq_len(nsim), \(i) {
        boot_mod <- bootstrap_model(model, ids)
        fit <- fit_nhmm(
          boot_mod$model, init, init_sd = 0, restarts = 0, lambda = lambda, 
          method = method, bound = bound, control = control,
          control_restart = list(), control_mstep = control_mstep
        )
        if (fit$estimation_results$return_code >= 0) {
          fit$gammas <- permute_states(fit$gammas, gammas_mle)
        } else {
          fit$gammas <- NULL
        }
        p()
        list(gammas = fit$gammas, idx = boot_mod$idx)
      }, future.seed = TRUE
    )
    idx <- do.call(cbind, lapply(out, "[[", "idx"))
    out <- lapply(out, "[[", "gammas")
  } else {
    S <- model$n_states
    formula_pi <- model$initial_formula
    formula_A <- model$transition_formula
    formula_B <- model$emission_formula
    time_var <- model$time_variable
    id_var <- model$id_variable
    out <- future.apply::future_lapply(
      seq_len(nsim), \(i) {
        mod <- simulate_nhmm(
          S, formula_B, formula_pi, formula_A,  model$data, id_var, time_var,
          init, init_sd = 0
        )$model
        fit <- fit_nhmm(
          mod, init, init_sd = 0, restarts = 0, lambda = lambda, 
          method = method, bound = bound, control = control,
          control_restart = list(), control_mstep = control_mstep
        )
        if (fit$estimation_results$return_code >= 0) {
          fit$gammas <- permute_states(fit$gammas, gammas_mle)
        } else {
          fit$gammas <- NULL
        }
        p()
        fit$gammas
      }, future.seed = TRUE
    )
  }
  boot <- list(
    gamma_pi = lapply(out, "[[", "gamma_pi"), 
    gamma_A = lapply(out, "[[", "gamma_A"), 
    gamma_B = lapply(out, "[[", "gamma_B")
  )
  boot <- lapply(boot,  function(x) x[lengths(x) > 0])
  if (length(boot[[1]]) < nsim) {
    warning_(
      paste0(
        "Estimation in some of the bootstrap samples failed. ",
        "Returning samples from {length(boot[[1]])} successes out of {nsim} ",
        "bootstrap samples."
      )
    )
  }
  if (type == "nonparametric") {
    boot$idx <- idx
  } else {
    boot$idx <- matrix(seq_len(model$n_sequences), model$n_sequences, nsim)
  }
  if (append && !is.null(model$boot)) {
    model$boot$gamma_pi <- c(model$boot$gamma_pi, boot$gamma_pi)
    model$boot$gamma_A <- c(model$boot$gamma_A, boot$gamma_A)
    model$boot$gamma_B <- c(model$boot$gamma_B, boot$gamma_B)
    model$boot$idx <- cbind(model$boot$idx, boot$idx)
  } else {
    model$boot <- boot
  }
  
  model
}
#' @rdname bootstrap
#' @export
bootstrap_coefs.mnhmm <- function(model, nsim, 
                                  type = c("nonparametric", "parametric"),
                                  append = FALSE, ...) {
  type <- try(match.arg(type), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val nonparametric} or 
    {.val parametric}."
  )
  stopifnot_(
    !missing(nsim) && checkmate::test_int(x = nsim, lower = 0L), 
    "Argument {.arg nsim} must be a single positive integer."
  )
  stopifnot_(
    checkmate::test_logical(x = append), 
    "Argument {.arg append} must be a single logical value."
  )
  init <- model$etas
  gammas_mle <- model$gammas
  D <- model$n_clusters
  pcp_mle <- matrix(
    posterior_cluster_probabilities(model)$probability, ncol = D, byrow = TRUE
  )
  lambda <- model$estimation_results$lambda
  bound <- model$estimation_results$bound
  method <- model$estimation_results$method
  p <- progressr::progressor(along = seq_len(nsim))
  original_options <- options(future.globals.maxSize = Inf)
  on.exit(options(original_options))
  control <- model$controls$control
  control$print_level <- 0
  control_mstep <- model$controls$mstep
  control_mstep$print_level <- 0
  if (type == "nonparametric") {
    ids <- unique(model$data[[model$id_variable]])
    out <- future.apply::future_lapply(
      seq_len(nsim), \(i) {
        boot_mod <- bootstrap_model(model, ids)
        fit <- fit_mnhmm(
          boot_mod$model, init, init_sd = 0, restarts = 0, lambda = lambda, 
          method = method, bound = bound, control = control,
          control_restart = list(), control_mstep = control_mstep
        )
        if (fit$estimation_results$return_code >= 0) {
          fit <- permute_clusters(fit, model, pcp_mle)
          for (j in seq_len(D)) {
            out <- permute_states(
              lapply(fit$gammas[c("gamma_pi", "gamma_A", "gamma_B")], "[[", j), 
              lapply(gammas_mle[c("gamma_pi", "gamma_A", "gamma_B")], "[[", j)
            )
            fit$gammas$gamma_pi[[j]] <- out$gamma_pi
            fit$gammas$gamma_A[[j]] <- out$gamma_A
            fit$gammas$gamma_B[[j]] <- out$gamma_B
          }
        } else {
          fit$gammas <- NULL
        }
        p()
        list(gammas = fit$gammas, idx = boot_mod$idx)
      }, future.seed = TRUE
    )
    idx <- do.call(cbind, lapply(out, "[[", "idx"))
    out <- lapply(out, "[[", "gammas")
  } else {
    S <- model$n_states
    formula_pi <- model$initial_formula
    formula_A <- model$transition_formula
    formula_B <- model$emission_formula
    formula_omega <- model$cluster_formula
    time_var <- model$time_variable
    id_var <- model$id_variable
    out <- future.apply::future_lapply(
      seq_len(nsim), \(i) {
        mod <- simulate_mnhmm(
          S, D, formula_B, formula_pi, formula_A, formula_omega, model$data, 
          id_var, time_var, init, init_sd = 0
        )$model
        fit <- fit_mnhmm(
          mod, init, init_sd = 0, restarts = 0, lambda = lambda, 
          method = method, bound = bound, control = control,
          control_restart = list(), control_mstep = control_mstep
        )
        if (fit$estimation_results$return_code >= 0) {
          fit <- permute_clusters(fit, model, pcp_mle)
          for (j in seq_len(D)) {
            out <- permute_states(
              lapply(fit$gammas[c("gamma_pi", "gamma_A", "gamma_B")], "[[", j), 
              lapply(gammas_mle[c("gamma_pi", "gamma_A", "gamma_B")], "[[", j)
            )
            fit$gammas$gamma_pi[[j]] <- out$gamma_pi
            fit$gammas$gamma_A[[j]] <- out$gamma_A
            fit$gammas$gamma_B[[j]] <- out$gamma_B
          }
        } else {
          fit$gammas <- NULL
        }
        p()
        fit$gammas
      }, future.seed = TRUE
    )
    
  }
  boot <- list(
    gamma_pi = lapply(out, "[[", "gamma_pi"), 
    gamma_A = lapply(out, "[[", "gamma_A"), 
    gamma_B = lapply(out, "[[", "gamma_B"),
    gamma_omega = lapply(out, "[[", "gamma_omega")
  )
  boot <- lapply(boot,  function(x) x[lengths(x) > 0])
  if (length(boot[[1]]) < nsim) {
    warning_(
      paste0(
        "Estimation in some of the bootstrap samples failed. ",
        "Returning samples from {length(boot[[1]])} successes out of {nsim} ",
        "bootstrap samples."
      )
    )
  }
  if (type == "nonparametric") {
    boot$idx <- idx
  } else {
    boot$idx <- matrix(seq_len(model$n_sequences), model$n_sequences, nsim)
  }
  if (append && !is.null(model$boot)) {
    model$boot$gamma_pi <- c(model$boot$gamma_pi, boot$gamma_pi)
    model$boot$gamma_A <- c(model$boot$gamma_A, boot$gamma_A)
    model$boot$gamma_B <- c(model$boot$gamma_B, boot$gamma_B)
    model$boot$gamma_omega <- c(model$boot$gamma_omega, boot$gamma_omega)
    model$boot$idx <- cbind(model$boot$idx, boot$idx)
  } else {
    model$boot <- boot
  }
  model
}
