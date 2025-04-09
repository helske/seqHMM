em_mnhmm <- function(model, inits, init_sd, restarts, lambda, 
                     bound, control, control_restart, control_mstep, 
                     save_all_solutions) {
  M <- model$n_symbols
  S <- model$n_states
  T_ <- model$length_of_sequences
  C <- model$n_channels
  D <- model$n_clusters
  np_pi <- attr(model, "np_pi")
  np_A <- attr(model, "np_A")
  np_B <- attr(model, "np_B")
  np_omega <- attr(model, "np_omega")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  X_omega <- model$X_omega
  icpt_only_pi <- io(X_pi)
  icpt_only_A <- io(X_A)
  icpt_only_B <- io(X_B)
  icpt_only_omega <- io(X_omega)
  iv_A <- iv(X_A)
  iv_B <- iv(X_B)
  tv_A <- tv(X_A)
  tv_B <- tv(X_B)
  K_omega <- nrow(X_omega)
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- vapply(model$X_B, \(x) nrow(x), 1L)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  obs <- create_obsArray(model)
  all_solutions <- NULL
  if (restarts > 0L) {
    .fun <- function(i, base_init, u) {
      init <- base_init
      z <- matrix(unlist(init$eta_pi) + u[seq_len(np_pi)], ncol = D)
      init$eta_pi <- lapply(seq_len(D), \(d) array(z[, d], dim = dim(init$eta_pi[[1]])))
      z <- matrix(unlist(init$eta_A) + u[np_pi + seq_len(np_A)], ncol = D)
      init$eta_A <- lapply(seq_len(D), \(d) array(z[, d], dim = dim(init$eta_A[[1]])))
      z <- array(unlist(init$eta_B) + u[np_pi + np_A + seq_len(np_B)], dim = c(np_B / (C * D), C, D))
      init$eta_B <- lapply(seq_len(D), \(d) lapply(seq_len(C), \(c) array(z[, c, d], dim = dim(init$eta_B[[1]][[c]]))))
      init$eta_omega <- init$eta_omega + u[np_pi + np_A + np_B + seq_len(np_omega)]
      fit <- EM_LBFGS_mnhmm(
        obs, Ti, M, X_pi, X_A, X_B, X_omega, 
        icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
        iv_A, iv_B, tv_A, tv_B, 
        init$eta_pi, init$eta_A, init$eta_B, init$eta_omega, lambda,
        control$maxeval, control$ftol_abs, control$ftol_rel,
        control$xtol_abs, control$xtol_rel, control$print_level, 
        control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, bound
      )
      p()
      fit
    }
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    base_init <- create_initial_values(inits, model, init_sd = 0)
    u <- t(
      stats::qnorm(
        lhs::maximinLHS(restarts, length(unlist(base_init))), sd = init_sd
      )
    )
    out <- future.apply::future_lapply(
      seq_len(restarts), \(i) .fun(i, base_init, u[, i]), future.seed = TRUE
    )
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
  out <- EM_LBFGS_mnhmm(
    obs, Ti, M, X_pi, X_A, X_B, X_omega, 
    icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
    iv_A, iv_B, tv_A, tv_B, 
    init$eta_pi, init$eta_A, init$eta_B, init$eta_omega, lambda,
    control$maxeval, control$ftol_abs, control$ftol_rel,
    control$xtol_abs, control$xtol_rel, control$print_level, 
    control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
    control_mstep$xtol_abs, control_mstep$xtol_rel, 
    control_mstep$print_level, bound
  )
  if (out$return_code < 0) {
    warning_(
      paste("Optimization terminated due to error:", return_msg(out$return_code))
    )
  }
  model$etas$pi <- drop(out$eta_pi)
  model$gammas$pi <- drop(eta_to_gamma_mat_field(model$etas$pi))
  model$etas$A <- drop(out$eta_A)
  model$gammas$A <- drop(eta_to_gamma_cube_field(model$etas$A))
  model$etas$B <- split(out$eta_B, seq_len(D))
  model$gammas$B <- split(eta_to_gamma_cube_2d_field(model$etas$B), seq_len(D))
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
