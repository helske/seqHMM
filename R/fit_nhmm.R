#' Estimate a Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_nhmm <- function(model, inits, init_sd, restarts, threads, verbose, penalize, penalty, ...) {
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  stopifnot_(
    checkmate::test_int(x = restarts, lower = 0L), 
    "Argument {.arg restarts} must be a single integer."
  )
  obs <- create_obsArray(model) + 1L
  if (model$n_channels == 1) {
    obs <- array(obs, dim(obs)[2:3])
  }
  K_i <- dim(model$X_initial)[2]
  K_s <- dim(model$X_transition)[3]
  K_o <- dim(model$X_emission)[3]
  M <- model$n_symbols
  S <- model$n_states
  if (identical(inits, "random")) {
    inits <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL)
  } else {
    if (is.null(inits$initial_probs)) inits$initial_probs <- NULL
    if (is.null(inits$transition_probs)) inits$transition_probs <- NULL
    if (is.null(inits$emission_probs)) inits$emission_probs <- NULL
  }
  model_code <- stanmodels[[attr(model, "type")]]
  if (restarts > 0L) {
    if (threads > 1L) {
      plan(multisession, workers = threads)
    } else {
      plan(sequential)
    }
    dots <- list(...)
    if (is.null(dots$hessian)) dots$hessian <- FALSE
    if (is.null(dots$tol_obj)) dots$tol_obj <- 1e-6
    if (is.null(dots$tol_rel_obj)) dots$tol_rel_obj <- 1e6
    if (is.null(dots$tol_grad)) dots$tol_grad <-  1e-4
    if (is.null(dots$tol_rel_grad)) dots$tol_rel_grad <- 1e9
    if (is.null(dots$tol_param)) dots$tol_param <- 1e-4
    if (is.null(dots$check_data)) dots$check_data <- FALSE
    out <- future_lapply(seq_len(restarts), function(i) {
      init <- create_initial_values(
        inits, S, M, init_sd, K_i, K_s, K_o
      )
      do.call(
        optimizing, 
        c(list(
          model_code, init = init,
          data = list(
            penalty = penalty,
            penalize = as.integer(penalize),
            N = model$n_sequences,
            T = model$sequence_lengths,
            max_T = model$length_of_sequences, 
            max_M = max(model$n_symbols),
            M = M,
            S = S,
            C = model$n_channels,
            K_i = K_i,
            K_s = K_s,
            K_o = K_o,
            X_i = model$X_initial,
            X_s = model$X_transition,
            X_o = model$X_emission,
            obs = obs,
            ids = seq_len(model$n_sequences),
            N_sample = model$n_sequences
          ),
          as_vector = FALSE,
          verbose = FALSE,
          save_iterations = FALSE
        ), dots)
      )[c("par", "value", "return_code")]
    },
    future.seed = TRUE)
    logliks <- unlist(lapply(out, "[[", "value"))
    return_codes <- unlist(lapply(out, "[[", "return_code"))
    successful <- which(return_codes == 0)
    optimum <- successful[which.max(logliks[successful])]
    init <- as.list(out[[optimum]]$par)
  } else {
    init <- create_initial_values(
      inits, S, M, init_sd, K_i, K_s, K_o
    )
  }
  dots <- list(...)
  if (is.null(dots$hessian)) dots$hessian <- TRUE
  if (is.null(dots$check_data)) dots$check_data <- FALSE
  if (is.null(dots$tol_obj)) dots$tol_obj <- 1e-12
  if (is.null(dots$tol_rel_obj)) dots$tol_rel_obj <- 1e4
  if (is.null(dots$tol_grad)) dots$tol_grad <-  1e-8
  if (is.null(dots$tol_rel_grad)) dots$tol_rel_grad <- 1e4
  if (is.null(dots$tol_param)) dots$tol_param <- 1e-8
  out <-  do.call(
    optimizing, 
    c(list(
      model_code, 
      data = list(
        penalty = penalty, 
        penalize = as.integer(penalize),
        N = model$n_sequences,
        T = model$sequence_lengths,
        max_T = model$length_of_sequences, 
        max_M = max(model$n_symbols),
        M = M,
        S = S,
        C = model$n_channels,
        K_i = K_i,
        K_s = K_s,
        K_o = K_o,
        X_i = model$X_initial,
        X_s = model$X_transition,
        X_o = model$X_emission,
        obs = obs,
        ids = seq_len(model$n_sequences),
        N_sample = model$n_sequences
      ), 
      as_vector = FALSE,
      init = init,
      verbose = verbose,
      save_iterations = FALSE
    ), dots)
  )[c("par", "value", "return_code", "hessian")]
  
  model$coefficients <- out$par[c("beta_i_raw", "beta_s_raw", "beta_o_raw")]
  model$stan_model <- model_code@model_code
  model$estimation_results <- list(
    hessian = out$hessian,
    penalized_loglik_N = out$par[["ploglik_N"]],
    penalized_loglik = out$value, 
    loglik = out$par[["log_lik"]], 
    penalty = out$par[["prior"]],
    return_code = out$return_code,
    plogliks_of_restarts = if(restarts > 1L) logliks else NULL, 
    return_codes_of_restarts = if(restarts > 1L) return_codes else NULL
  )
  model
}
