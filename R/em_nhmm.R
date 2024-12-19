em_nhmm <- function(model, inits, init_sd, restarts, lambda, 
                    bound, control, control_restart, control_mstep, 
                    save_all_solutions) {
  M <- model$n_symbols
  S <- model$n_states
  T_ <- model$length_of_sequences
  C <- model$n_channels
  icpt_only_pi <- attr(model$X_pi, "icpt_only")
  icpt_only_A <- attr(model$X_A, "icpt_only")
  icpt_only_B <- attr(model$X_B, "icpt_only")
  iv_A <- attr(model$X_A, "iv")
  iv_B <- attr(model$X_B, "iv")
  tv_A <- attr(model$X_A, "tv")
  tv_B <- attr(model$X_B, "tv")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  obs <- create_obsArray(model)
  if (C == 1) {
    obs <- array(obs, dim(obs)[2:3])
  }
  all_solutions <- NULL
  if (restarts > 0L) {
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    out <- future.apply::future_lapply(seq_len(restarts), function(i) {
      init <- create_initial_values(inits, model, init_sd)
      if (C == 1) {
        fit <- EM_LBFGS_nhmm_singlechannel(
          init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, obs,
          Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
          n_obs, control_restart$maxeval,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, lambda, bound)
      } else {
        fit <- EM_LBFGS_nhmm_multichannel(
          init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, obs,
          Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
          n_obs, control_restart$maxeval,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, lambda, bound)
      }
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
    init <- optimum[c("eta_pi", "eta_A", "eta_B")]
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
  }
  if (C == 1) {
    out <- EM_LBFGS_nhmm_singlechannel(
      init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, obs,
      Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
      n_obs, control$maxeval,
      control$ftol_abs, control$ftol_rel,
      control$xtol_abs, control$xtol_rel, 
      control$print_level, control_mstep$maxeval,
      control_mstep$ftol_abs, control_mstep$ftol_rel,
      control_mstep$xtol_abs, control_mstep$xtol_rel, 
      control_mstep$print_level, lambda, bound)
  } else {
    out <- EM_LBFGS_nhmm_multichannel(
      init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, obs,
      Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
      n_obs, control$maxeval,
      control$ftol_abs, control$ftol_rel,
      control$xtol_abs, control$xtol_rel, 
      control$print_level, control_mstep$maxeval,
      control_mstep$ftol_abs, control_mstep$ftol_rel,
      control_mstep$xtol_abs, control_mstep$xtol_rel, 
      control_mstep$print_level, lambda, bound)
  }
  if (out$return_code < 0) {
    warning_(
      paste("Optimization terminated due to error:", return_msg(out$return_code))
    )
  }
  model$etas$pi[] <- out$eta_pi
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$etas$A <- out$eta_A
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  if (C == 1L) {
    model$etas$B[] <- out$eta_B
    model$gammas$B <- eta_to_gamma_cube(model$etas$B)
  } else {
    for (i in seq_len(C)) model$etas$B[[i]] <- out$eta_B[[i]]
    model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
  }
  
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
