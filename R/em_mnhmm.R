em_mnhmm <- function(model, inits, init_sd, restarts, lambda, 
                     bound, control, control_restart, control_mstep, 
                     save_all_solutions) {
  M <- model$n_symbols
  S <- model$n_states
  T_ <- model$length_of_sequences
  C <- model$n_channels
  D <- model$n_clusters
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  n_d <- attr(model, "np_omega")
  icpt_only_omega <- attr(model$X_omega, "icpt_only")
  icpt_only_pi <- attr(model$X_pi, "icpt_only")
  icpt_only_A <- attr(model$X_A, "icpt_only")
  icpt_only_B <- attr(model$X_B, "icpt_only")
  iv_A <- attr(model$X_A, "iv")
  iv_B <- attr(model$X_B, "iv")
  tv_A <- attr(model$X_A, "tv")
  tv_B <- attr(model$X_B, "tv")
  X_omega <- model$X_omega
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  K_omega <- nrow(X_omega)
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  obs <- create_obsArray(model)
  all_solutions <- NULL
  if (restarts > 0L) {
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    base_init <- create_initial_values(inits, model, init_sd = 0)
    u <- t(stats::qnorm(
      lhs::maximinLHS(restarts, length(unlist(base_init))), 
      sd = init_sd))
    out <- future.apply::future_lapply(seq_len(restarts), function(i) {
      init <- base_init
      init$eta_pi[] <- init$eta_pi[] + u[seq_len(n_i), i]
      init$eta_A[] <- init$eta_A[] + u[n_i + seq_len(n_s), i]
      init$eta_B[] <- init$eta_B[] + u[n_i + n_s + seq_len(n_o), i]
      init$eta_omega[] <- init$eta_omega[] + u[n_i + n_s + n_o + seq_len(n_d), i]
      if (C == 1) {
        fit <- EM_LBFGS_mnhmm_singlechannel(
          init$eta_omega, model$X_omega, init$eta_pi, model$X_pi, init$eta_A, model$X_A, 
          init$eta_B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
          icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
          n_obs, control_restart$maxeval,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, lambda, bound)
      } else {
        eta_B <- unlist(init$eta_B, recursive = FALSE)
        fit <- EM_LBFGS_mnhmm_multichannel(
          init$eta_omega, model$X_omega, init$eta_pi, model$X_pi, init$eta_A, model$X_A, 
          eta_B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
          icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
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
      optimum <- out[[1]]
    } else {
      logliks <- unlist(lapply(out, "[[", "logLik"))
      optimum <- out[[which.max(logliks)]]
    }
    init <- optimum[c("eta_omega", "eta_pi", "eta_A", "eta_B")]
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
  }
  if (C == 1) {
    out <- EM_LBFGS_mnhmm_singlechannel(
      init$eta_omega, model$X_omega, init$eta_pi, model$X_pi, init$eta_A, model$X_A, 
      init$eta_B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
      icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
      n_obs, control$maxeval,
      control$ftol_abs, control$ftol_rel,
      control$xtol_abs, control$xtol_rel, 
      control$print_level, control_mstep$maxeval,
      control_mstep$ftol_abs, control_mstep$ftol_rel,
      control_mstep$xtol_abs, control_mstep$xtol_rel, 
      control_mstep$print_level, lambda, bound)
  } else {
    eta_B <- unlist(init$eta_B, recursive = FALSE)
    out <- EM_LBFGS_mnhmm_multichannel(
      init$eta_omega, model$X_omega, init$eta_pi, model$X_pi, init$eta_A, model$X_A, 
      eta_B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
      icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
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
  model$etas$pi <- c(out$eta_pi)
  model$gammas$pi <- c(eta_to_gamma_mat_field(model$etas$pi))
  model$etas$A <- c(out$eta_A)
  model$gammas$A <- c(eta_to_gamma_cube_field(model$etas$A))
  if (C == 1L) {
    model$etas$B <- c(out$eta_B)
    model$gammas$B <- c(eta_to_gamma_cube_field(model$etas$B))
  } else {
    l <- lengths(model$etas$B)
    model$etas$B <- unname(split(c(out$eta_B), rep(seq_along(l), l)))
    gamma_B <- c(eta_to_gamma_cube_field(unlist(model$etas$B, recursive = FALSE)))
    model$gammas$B <- unname(split(gamma_B, rep(seq_along(l), l)))
  }
  model$etas$omega <- out$eta_omega
  model$gammas$omega <- eta_to_gamma_mat(model$etas$omega)
  
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
