em_dnm_mnhmm <- function(model, inits, init_sd, restarts, lambda, 
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
  K_omega <- K(X_omega)
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  Ti <- model$sequence_lengths
  need_grad <- grepl("NLOPT_LD_", control$algorithm)
  obs <- create_obs(model)
  all_solutions <- NULL
  objectivef <- make_objective_mnhmm(
    model, lambda, grepl("NLOPT_LD_", control$algorithm)
  )
  use_fanhmm <- inherits(model, "fanhmm") && !identical(model$prior_obs, 0L)
  if (restarts > 0L) {
    .fun <- function(base_init, u) {
      pars <- base_init + u
      eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
      eta_A <- create_eta_A_mnhmm(
        pars[np_pi + seq_len(np_A)], 
        S, K_A, D
      )
      eta_B <- create_eta_B_mnhmm(
        pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
      )
      eta_omega <- create_eta_omega_mnhmm(
        pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
      )
      if (use_fanhmm) {
        fit <- Rcpp_EM_LBFGS_mfanhmm(
          obs, Ti, M, X_pi, X_A, X_B, X_omega, 
          icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
          iv_A, iv_B, tv_A, tv_B, 
          eta_pi, eta_A, eta_B, eta_omega, 
          model$prior_obs, model$W_X_B, lambda,
          control_restart$maxeval_em_dnm, control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, control_restart$print_level, 
          control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, bound, control_mstep$tolg
        )
      } else {
        fit <- Rcpp_EM_LBFGS_mnhmm(
          obs, Ti, M, X_pi, X_A, X_B, X_omega, 
          icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
          iv_A, iv_B, tv_A, tv_B, 
          eta_pi, eta_A, eta_B, eta_omega, lambda,
          control_restart$maxeval_em_dnm, control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, control_restart$print_level, 
          control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, bound, control_mstep$tolg
        )
      }
      em_return_code <- fit$return_code
      if (em_return_code >= 0) {
        init <- unlist(fit[c("eta_pi", "eta_A", "eta_B", "eta_omega")])
        fit <- nloptr(
          x0 = init, eval_f = objectivef, lb = -rep(bound, length(init)), 
          ub = rep(bound, length(init)), opts = control_restart
        )
        p()
        fit
      } else {
        p()
        list(status = -1000 + em_return_code, objective = NaN)
      }
    }
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    base_init <- unlist(create_initial_values(inits, model, init_sd = 0))
    u <- t(
      stats::qnorm(
        lhs::maximinLHS(restarts, length(unlist(base_init))), 
        sd = init_sd
      )
    )
    fit <- future.apply::future_lapply(
      seq_len(restarts), \(i) .fun(base_init, u[, i]), future.seed = TRUE
    )
    logliks <- -unlist(lapply(fit, "[[", "objective")) * nobs(model)
    return_codes <- unlist(lapply(fit, "[[", "status"))
    successful <- which(return_codes > 0)
    if (length(successful) == 0) {
      warning_(
        c("All restarts terminated due to error.",
          "Error of first restart: ", return_msg(return_codes[1]),
          "Running DNM using initial values for EM.")
      )
      init <- create_initial_values(inits, model, init_sd)
      em_return_code <- return_codes[1] + 1000
    } else {
      em_return_code <- 0 # generic success
      optimum <- successful[which.max(logliks[successful])]
      pars <- fit[[optimum]]$solution
      init <- list()
      init$eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
      init$eta_A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
      init$eta_B <- create_eta_B_mnhmm(
        pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
      )
      init$eta_omega <- create_eta_omega_mnhmm(
        pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
      )
    }
    if (save_all_solutions) {
      all_solutions <- fit
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
    if (use_fanhmm) {
      fit <- Rcpp_EM_LBFGS_mfanhmm(
        obs, Ti, M, X_pi, X_A, X_B, X_omega, 
        icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
        iv_A, iv_B, tv_A, tv_B, 
        init$eta_pi, init$eta_A, init$eta_B, init$eta_omega, 
        model$prior_obs, model$W_X_B, lambda,
        control$maxeval_em_dnm, control$ftol_abs, control$ftol_rel,
        control$xtol_abs, control$xtol_rel, control$print_level, 
        control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, bound, control_mstep$tolg
      )
    } else {
      fit <- Rcpp_EM_LBFGS_mnhmm(
        obs, Ti, M, X_pi, X_A, X_B, X_omega, 
        icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
        iv_A, iv_B, tv_A, tv_B, 
        init$eta_pi, init$eta_A, init$eta_B, init$eta_omega, lambda,
        control$maxeval_em_dnm, control$ftol_abs, control$ftol_rel,
        control$xtol_abs, control$xtol_rel, control$print_level, 
        control_mstep$maxeval, control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, bound, control_mstep$tolg
      )
    }
    em_return_code <- fit$return_code
    if (em_return_code >= 0) {
      init <- fit[c("eta_pi", "eta_A", "eta_B", "eta_omega")]
    } else {
      warning_(
        paste("EM-step terminated due to error:", return_msg(em_return_code),
              "Running DNM using initial values for EM.")
      )
      init <- create_initial_values(inits, model, init_sd)
    }
  }
  fit <- nloptr(
    x0 = unlist(init), eval_f = objectivef, 
    lb = -rep(bound, length(unlist(init))), 
    ub = rep(bound, length(unlist(init))), opts = control
  )
  if (fit$status < 0) {
    warning_(
      paste("Optimization terminated due to error:", return_msg(fit$status))
    )
    loglik <- NaN
  } else {
    loglik <- -fit$objective * nobs(model)
  }
  
  pars <- fit$solution
  model$etas$eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
  model$gammas$gamma_pi <- drop(eta_to_gamma_mat_field(model$etas$eta_pi))
  model$etas$eta_A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
  model$gammas$gamma_A <- drop(eta_to_gamma_cube_field(model$etas$eta_A))
  model$etas$eta_B <- create_eta_B_mnhmm(
    pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
  )
  model$gammas$gamma_B <- split(
    eta_to_gamma_cube_2d_field(model$etas$eta_B), seq_len(D)
  )
  model$etas$eta_omega <- create_eta_omega_mnhmm(
    pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
  )
  model$gammas$gamma_omega <- eta_to_gamma_mat(model$etas$eta_omega)
  model$estimation_results <- list(
    loglik = loglik, 
    return_code = fit$status,
    message = fit$message,
    iterations = fit$iterations,
    logliks_of_restarts = if(restarts > 0L) logliks else NULL, 
    return_codes_of_restarts = if(restarts > 0L) return_codes else NULL,
    all_solutions = all_solutions,
    lambda = lambda,
    bound = bound,
    method = "EM-DNM",
    algorithm = control$algorithm,
    EM_return_code = em_return_code
  )
  model
}
