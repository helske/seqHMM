#' Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_mnhmm <- function(model, inits, init_sd, restarts, threads, ...) {
  check_positive_integer(threads)
  obs <- create_obsArray(model) + 1L
  if (model$n_channels == 1) {
    obs <- obs[1, ,]
  }
  K_i <- dim(model$X_initial)[2]
  K_s <- dim(model$X_transition)[3]
  K_o <- dim(model$X_emission)[3]
  K_d <- dim(model$X_cluster)[2]
  D <- model$n_clusters
  M <- model$n_symbols
  S <- model$n_states
  if (identical(inits, "random")) {
    inits <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL,
      cluster_probs = NULL)
  } else {
    if (is.null(inits$initial_probs)) inits$initial_probs <- NULL
    if (is.null(inits$transition_probs)) inits$transition_probs <- NULL
    if (is.null(inits$emission_probs)) inits$emission_probs <- NULL
    if (is.null(inits$cluster_probs)) inits$cluster_probs <- NULL
  }
  if (restarts > 1L) {
    if (threads > 1L) {
      plan(multisession, workers = threads)
    } else {
      plan(sequential)
    }
    out <- future_lapply(seq_len(restarts), function(i) {
      init <- create_initial_values(
        inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
      )
      optimizing(
        stanmodels[[attr(model, "type")]], init = init,
        data = list(
          N = model$n_sequences,
          T = model$length_of_sequences, 
          max_M = max(model$n_symbols),
          M = M,
          S = S,
          C = model$n_channels,
          D = D,
          K_i = K_i,
          K_s = K_s,
          K_o = K_o,
          K_d = K_d,
          X_i = model$X_initial,
          X_s = model$X_transition,
          X_o = model$X_emission,
          X_d = model$X_cluster,
          obs = obs),
        check_data = FALSE,
        as_vector = FALSE,
        tol_obj = 1e-6, tol_rel_obj = 1e6,
        tol_grad = 1e-4, tol_rel_grad = 1e9,
        tol_param = 1e-4,
        ...)[c("par", "value", "return_code")]
    },
    future.seed = TRUE)
    logliks <- unlist(lapply(out, "[[", "value"))
    return_codes <- unlist(lapply(out, "[[", "return_code"))
    successful <- which(return_codes == 0)
    optimum <- successful[which.max(logliks[successful])]
    init <- as.list(out[[optimum]]$par)
  } else {
    init <- create_initial_values(
      inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
    )
  }
  out <- optimizing(
    stanmodels[[attr(model, "type")]], 
    data = list(
      N = model$n_sequences,
      T = model$length_of_sequences, 
      max_M = max(model$n_symbols),
      M = M,
      S = S,
      C = model$n_channels,
      D = D,
      K_i = K_i,
      K_s = K_s,
      K_o = K_o,
      K_d = K_d,
      X_i = model$X_initial,
      X_s = model$X_transition,
      X_o = model$X_emission,
      X_d = model$X_cluster,
      obs = obs), 
    check_data = FALSE,
    as_vector = FALSE,
    hessian = TRUE,
    init = init,
    ...)[c("par", "value", "return_code", "hessian")]
  
  model$estimation_results <- list(
    parameters = out$par[lengths(out$par) > 0], 
    hessian = out$hessian,
    logLik = out$value, 
    return_code = out$return_code,
    logLiks_of_restarts = if(restarts > 1L) logliks else NULL, 
    return_codes_of_restarts = if(restarts > 1L) return_codes else NULL
  )
  model
}
