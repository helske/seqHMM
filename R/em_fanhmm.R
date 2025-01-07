em_fanhmm <- function(model, inits, init_sd, restarts, lambda, 
                      bound, control, control_restart, control_mstep, 
                      save_all_solutions) {
  M <- model$n_symbols
  S <- model$n_states
  T_ <- model$length_of_sequences
  #  C <- model$n_channels
  icpt_only_pi <- attr(model$X_pi, "icpt_only")
  icpt_only_A <- attr(model$X_A, "icpt_only")
  iv_A <- attr(model$X_A, "iv")
  iv_B <- attr(model$X_B, "iv")
  tv_A <- attr(model$X_A, "tv")
  tv_B <- attr(model$X_B, "tv")
  if (!is.null(model$feedback_formula)) {
    icpt_only_A <- FALSE
    iv_A <- TRUE
    tv_A <- TRUE
  }
  icpt_only_B <- attr(model$X_B, "icpt_only")
  if (!is.null(model$autoregression_formula)) {
    icpt_only_B <- FALSE
    iv_B <- TRUE
    tv_B <- TRUE
  }
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B)
  W_A <- model$W_A
  W_B <- model$W_B
  L_A <- nrow(W_A)
  L_B <- nrow(W_B)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  obs <- create_obsArray(model)
  obs <- array(obs, dim(obs)[2:3])
  all_solutions <- NULL
  if (restarts > 0L) {
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    out <- future.apply::future_lapply(seq_len(restarts), function(i) {
      init <- c(
        create_initial_values(inits, model, init_sd),
        rho_A = list(create_rho_A_inits(inits$rho_A, S, M, L_A, init_sd)), 
        rho_B = list(create_rho_B_inits(inits$rho_B, S, M, L_B, init_sd))
      )
      fit <- EM_LBFGS_fanhmm_singlechannel(
        init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, 
        init$rho_A, model$W_A, init$rho_B, model$W_B, model$obs_0, obs,
        Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
        n_obs, control_restart$maxeval,
        control_restart$ftol_abs, control_restart$ftol_rel,
        control_restart$xtol_abs, control_restart$xtol_rel, 
        control_restart$print_level, control_mstep$maxeval,
        control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, lambda, bound)
      p()
      fit
    },
    future.seed = TRUE)
    return_codes <- unlist(lapply(out, "[[", "return_code"))
    if (all(return_codes < 0)) {
      warning_(
        c("All restarts terminated due to error.",
          "Error of first restart: ", return_msg(return_codes[1]))
      )
    }
    logliks <- unlist(lapply(out, "[[", "logLik"))
    optimum <- out[[which.max(logliks)]]
    init <- optimum[c("eta_pi", "eta_A", "eta_B", "rho_A", "rho_B")]
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- c(
      create_initial_values(inits, model, init_sd),
      rho_A = list(create_rho_A_inits(inits$rho_A, S, M, L_A, init_sd)), 
      rho_B = list(create_rho_B_inits(inits$rho_B, S, M, L_B, init_sd))
    )
  }
  out <- EM_LBFGS_fanhmm_singlechannel(
    init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, 
    init$rho_A, model$W_A, init$rho_B, model$W_B, model$obs_0, obs,
    Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
    n_obs, control$maxeval,
    control$ftol_abs, control$ftol_rel,
    control$xtol_abs, control$xtol_rel, 
    control$print_level, control_mstep$maxeval,
    control_mstep$ftol_abs, control_mstep$ftol_rel,
    control_mstep$xtol_abs, control_mstep$xtol_rel, 
    control_mstep$print_level, lambda, bound)
  if (out$return_code < 0) {
    warning_(
      paste("Optimization terminated due to error:", return_msg(out$return_code))
    )
  }
  model$etas$pi[] <- out$eta_pi
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$etas$A <- out$eta_A
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  model$etas$B[] <- out$eta_B
  model$gammas$B <- eta_to_gamma_cube(model$etas$B)
  model$rhos$A <- out$rho_A
  model$phis$A <- rho_to_phi_field(model$rhos$A)
  model$rhos$B <- out$rho_B
  model$phis$B <- rho_to_phi_field(model$rhos$B)
  model$estimation_results <- list(
    loglik = out$logLik,
    penalty = out$penalty_term,
    iterations = out$iterations,
    return_code = out$return_code,
    logliks_of_restarts = if(restarts > 0L) logliks else NULL, 
    return_codes_of_restarts = if(restarts > 0L) return_codes else NULL,
    all_solutions = all_solutions,
    f_rel_change = out$relative_f_change,
    f_abs_change = out$absolute_f_change,
    x_rel_change = out$relative_x_change,
    x_abs_change = out$absolute_x_change,
    lambda = lambda,
    bound = bound,
    method = "EM"
  )
  model
}
