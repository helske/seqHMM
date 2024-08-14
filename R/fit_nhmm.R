#' Estimate Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_nhmm <- function(model, restarts, threads, ...) {
  
  obs <- t(sapply(model$observations, as.integer))
  obs[obs > model$n_symbols] <- 0L # missing values to zero
  
  if (restarts > 1L) {
    if (threads > 1L) {
      plan(multisession, workers = threads)
    } else {
      plan(sequential)
    }
    out <- future_lapply(seq_len(restarts), function(i) {
      optimizing(
        stanmodels[[attr(model, "type")]], init = "random",
        data = list(
          N = model$n_sequences,
          T = model$length_of_sequences, 
          M = model$n_symbols,
          S = model$n_states,
          C = model$n_channels,
          K_i = dim(model$X_initial)[2],
          K_s = dim(model$X_transition)[3],
          K_o = dim(model$X_emission)[3],
          X_i = model$X_initial,
          X_s = model$X_transition,
          X_o = model$X_emission,
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
    init <- "random"
  }
  
  out <- optimizing(
    stanmodels[[attr(model, "type")]], 
    data = list(
      N = model$n_sequences,
      T = model$length_of_sequences, 
      M = model$n_symbols,
      S = model$n_states,
      C = model$n_channels,
      K_i = dim(model$X_initial)[2],
      K_s = dim(model$X_transition)[3],
      K_o = dim(model$X_emission)[3],
      X_i = model$X_initial,
      X_s = model$X_transition,
      X_o = model$X_emission,
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
