em_nhmm <- function(model, inits, init_sd, restarts, lambda, 
                    bound, control, control_restart, control_mstep, 
                    save_all_solutions) {
  M <- model$n_symbols
  S <- model$n_states
  T_ <- model$length_of_sequences
  C <- model$n_channels
  np_pi <- attr(model, "np_pi")
  np_A <- attr(model, "np_A")
  np_B <- attr(model, "np_B")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  icpt_only_pi <- io(X_pi)
  icpt_only_A <- io(X_A)
  icpt_only_B <- io(X_B)
  iv_A <- iv(X_A)
  iv_B <- iv(X_B)
  tv_A <- tv(X_A)
  tv_B <- tv(X_B)
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  Ti <- model$sequence_lengths
  obs <- create_obs(model)
  all_solutions <- NULL
  use_fanhmm <- inherits(model, "fanhmm") && !identical(model$prior_obs, 0L)
  if (restarts > 0L) {
    .fun <- function(base_init, u) {
      pars <- base_init + u
      eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
      eta_A <- create_eta_A_nhmm(
        pars[np_pi + seq_len(np_A)], 
        S, K_A
      )
      eta_B <- create_eta_B_nhmm(
        pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
      )
      if (use_fanhmm) {
        fit <- Rcpp_EM_LBFGS_fanhmm(
          obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B,
          iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, 
          model$prior_obs, model$W_X_B, lambda,
          control_restart$maxeval,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, bound, control_mstep$tolg
        )
      } else {
        fit <- Rcpp_EM_LBFGS_nhmm(
          obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B,
          iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, lambda,
          control_restart$maxeval,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, bound, control_mstep$tolg
        )
      }
      p()
      fit
    }
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    base_init <- unlist(create_initial_values(inits, model, init_sd = 0))
    u <- t(
      stats::qnorm(
        lhs::maximinLHS(restarts, length(unlist(base_init))), sd = init_sd
      )
    )
    fit <- future.apply::future_lapply(
      seq_len(restarts), \(i) .fun(base_init, u[, i]), future.seed = TRUE
    )
    return_codes <- unlist(lapply(fit, "[[", "return_code"))
    if (all(return_codes < 0)) {
      warning_(
        c("All restarts terminated due to error.",
          "Error of first restart: ", return_msg(return_codes[1]))
      )
      optimum <- fit[[1]]
    } else {
      logliks <- unlist(lapply(fit, "[[", "logLik"))
      optimum <- fit[[which.max(logliks)]]
    }
    init <- create_initial_values(
      optimum[c("eta_pi", "eta_A", "eta_B")], 
      model, 
      init_sd = 0
    )
    if (save_all_solutions) {
      all_solutions <- fit
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
  }
  if (use_fanhmm) {
    fit <- Rcpp_EM_LBFGS_fanhmm(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B,
      iv_A, iv_B, tv_A, tv_B, init$eta_pi, init$eta_A, init$eta_B, 
      model$prior_obs, model$W_X_B, lambda,
      control$maxeval, control$ftol_abs, control$ftol_rel,
      control$xtol_abs, control$xtol_rel, control$print_level, 
      control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
      control_mstep$xtol_abs, control_mstep$xtol_rel, 
      control_mstep$print_level, bound, control_mstep$tolg
    )
  } else {
    fit <- Rcpp_EM_LBFGS_nhmm(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B,
      iv_A, iv_B, tv_A, tv_B, init$eta_pi, init$eta_A, init$eta_B, lambda,
      control$maxeval, control$ftol_abs, control$ftol_rel,
      control$xtol_abs, control$xtol_rel, control$print_level, 
      control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
      control_mstep$xtol_abs, control_mstep$xtol_rel, 
      control_mstep$print_level, bound, control_mstep$tolg
    )
  }
  if (fit$return_code < 0) {
    warning_(
      paste("Optimization terminated due to error:", return_msg(fit$return_code))
    )
  }
  model$etas$eta_pi <- fit$eta_pi
  model$gammas$gamma_pi <- fit$gamma_pi
  model$etas$eta_A <- fit$eta_A
  model$gammas$gamma_A <- fit$gamma_A
  model$etas$eta_B <- drop(fit$eta_B)
  model$gammas$gamma_B <- drop(fit$gamma_B)
  
  model$estimation_results <- list(
    loglik = fit$logLik,
    iterations = fit$iterations,
    return_code = fit$return_code,
    logliks_of_restarts = if(restarts > 0L) logliks else NULL,
    return_codes_of_restarts = if(restarts > 0L) return_codes else NULL,
    all_solutions = all_solutions,
    f_rel_change = fit$relative_f_change,
    f_abs_change = fit$absolute_f_change,
    x_rel_change = fit$relative_x_change,
    x_abs_change = fit$absolute_x_change,
    lambda = lambda,
    bound = bound,
    method = "EM"
  )
  model
}
