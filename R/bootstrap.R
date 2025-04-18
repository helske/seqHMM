#' Bootstrap observations
#' @noRd
bootstrap_model <- function(model, ids) {
  id <- y <- NULL
  idx <- sort(sample.int(model$n_sequences, replace = TRUE))
  ids <- ids[idx]
  model$data <- model$data[list(ids), on = model$id_variable]
  new_ids <- rep(seq_len(model$n_sequences), each = model$length_of_sequences)
  model$data[, id := new_ids, env = list(id = model$id_variable)]
  model$sequence_lengths <- model$sequence_lengths[idx]
  model$X_pi[] <- model$X_pi[, idx]
  model$X_A[] <- model$X_A[, , idx]
  model$X_B[] <- model$X_B[, , idx]
  if (inherits(model, "mnhmm")) {
    model$X_omega[] <- model$X_omega[, idx]
  }
  lag_obs <- paste0("lag_", model$responses)
  ar <- any(
    lag_obs %in% attr(stats::terms(model$emission_formula), "term.labels")
  )
  attr(model, "nobs") <- sum(!is.na(
    model$data[, y, env = list(y = I(model$responses))]
  )) / model$n_channels - ar * model$n_sequences
  list(model = model, idx = idx)
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
  pcp <- matrix(
    posterior_cluster_probabilities(model)$probability, 
    ncol = model$n_clusters, byrow = TRUE
  )
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
#' @param nsim number of bootstrap samples.
#' @param type Either `"nonparametric"` (default) or `"parametric"`, to define 
#' whether nonparametric or parametric bootstrap should be used. The former samples 
#' sequences with replacement, whereas the latter simulates new datasets based 
#' on the model.
#' @param method Estimation method used in bootstrapping. Defaults to `"EM-DNM"`.
#' @param append If `TRUE`, in case the model already contains 
#' bootstrap samples, new samples are appended to `model$boot`. If `FALSE` 
#' (default), old samples are discarded.
#' @param ... Additional arguments to [estimate_nhmm()] or [estimate_mnhmm()].
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
  init <- stats::setNames(model$etas, c("eta_pi", "eta_A", "eta_B"))
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
      seq_len(nsim), function(i) {
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
    N <- model$n_sequences
    T_ <- model$sequence_lengths
    M <- model$n_symbols
    S <- model$n_states
    formula_pi <- model$initial_formula
    formula_A <- model$transition_formula
    formula_B <- model$emission_formula
    d <- copy(model$data)
    time_var <- model$time_variable
    id_var <- model$id_variable
    y <- model$responses
    d[, y := NULL, env = list(y = I(y))]
    out <- future.apply::future_lapply(
      seq_len(nsim), function(i) {
        mod <- simulate_nhmm(
          N, T_, M, S, formula_pi, formula_A, formula_B,
          d, id_var, time_var, init, 0, responses = y)$model
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
    gamma_pi = lapply(out, "[[", "pi"), 
    gamma_A = lapply(out, "[[", "A"), 
    gamma_B = lapply(out, "[[", "B")
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
  init <- stats::setNames(model$etas, c("eta_pi", "eta_A", "eta_B", "eta_omega"))
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
      seq_len(nsim), function(i) {
        boot_mod <- bootstrap_model(model, ids)
        fit <- fit_mnhmm(
          boot_mod$model, init, init_sd = 0, restarts = 0, lambda = lambda, 
          method = method, bound = bound, control = control,
          control_restart = list(), control_mstep = control_mstep
        )
        if (fit$estimation_results$return_code >= 0) {
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
    N <- model$n_sequences
    T_ <- model$sequence_lengths
    M <- model$n_symbols
    S <- model$n_states
    formula_pi <- model$initial_formula
    formula_A <- model$transition_formula
    formula_B <- model$emission_formula
    formula_omega <- model$cluster_formula
    d <- copy(model$data)
    time_var <- model$time_variable
    id_var <- model$id_variable
    y <- model$responses
    d[, y := NULL, env = list(y = I(y))]
    out <- future.apply::future_lapply(
      seq_len(nsim), function(i) {
        mod <- simulate_mnhmm(
          N, T_, M, S, D, formula_pi, formula_A, formula_B, formula_omega,
          d, id_var, time_var, init, 0, responses = y)$model
        fit <- fit_mnhmm(
          mod, init, init_sd = 0, restarts = 0, lambda = lambda, 
          method = method, bound = bound, control = control,
          control_restart = list(), control_mstep = control_mstep
        )
        if (fit$estimation_results$return_code >= 0) {
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
        } else {
          fit$gammas <- NULL
        }
        p()
        fit$gammas
      }, future.seed = TRUE
    )
    
  }
  boot <- list(
    gamma_pi = lapply(out, "[[", "pi"), 
    gamma_A = lapply(out, "[[", "A"), 
    gamma_B = lapply(out, "[[", "B"),
    gamma_omega = lapply(out, "[[", "omega")
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
