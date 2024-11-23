em_dnm_nhmm <- function(model, inits, init_sd, restarts, lambda, pseudocount, 
                          control, control_restart, control_mstep, 
                          save_all_solutions) {
  M <- model$n_symbols
  S <- model$n_states
  T_ <- model$length_of_sequences
  C <- model$n_channels
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
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
  need_grad <- grepl("NLOPT_LD_", control$algorithm)
  obs <- create_obsArray(model)
  if (C == 1) {
    obs <- array(obs, dim(obs)[2:3])
  }
  all_solutions <- NULL
  if (C == 1L) {
    if (need_grad) {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
        eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
        eta_B <- create_eta_B_nhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B
        )
        out <- log_objective_nhmm_singlechannel(
          eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, Ti,
          icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B
        )
        list(
          objective = - (out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs, 
          gradient = - (unlist(out[-1]) - lambda * pars) / n_obs
        )
      }
    } else {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
        eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
        eta_B <- create_eta_B_nhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B
        )
        out <- forward_nhmm_singlechannel(
          eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, Ti,
          icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B
        )
        - (sum(apply(out[, T_, ], 2, logSumExp)) - 0.5 * lambda * sum(pars^2)) / n_obs
      }
    }
  } else {
    if (need_grad) {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
        eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
        eta_B <- create_eta_multichannel_B_nhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B
        )
        out <- log_objective_nhmm_multichannel(
          eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, Ti,
          icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B
        )
        list(
          objective = - (out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs, 
          gradient = - (unlist(out[-1]) - lambda * pars) / n_obs
        )
      }
    } else {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
        eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
        eta_B <- create_eta_multichannel_B_nhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B
        )
        out <- forward_nhmm_multichannel(
          eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, Ti,
          icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B
        )
        - (sum(apply(out[, T_, ], 2, logSumExp)) - 0.5 * lambda * sum(pars^2)) / n_obs
      }
    }
  }
  
  start_time <- proc.time()
  if (restarts > 0L) {
    p <- progressr::progressor(along = seq_len(restarts))
    out <- future.apply::future_lapply(seq_len(restarts), function(i) {
      init <- create_initial_values(inits, model, init_sd)
      if (C == 1) {
        fit <- EM_LBFGS_nhmm_singlechannel(
          init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
          Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
          n_obs, control$maxeval_em_dnm,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, lambda, pseudocount)
      } else {
        fit <- EM_LBFGS_nhmm_multichannel(
          init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
          Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
          n_obs, control$maxeval_em_dnm,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, lambda, pseudocount)
      }
      if (fit$return_code == 0) {
        init <- unlist(
          create_initial_values(
            stats::setNames(
              fit[c("eta_pi", "eta_A", "eta_B")], c("pi", "A", "B")
            ), 
            model, 
            init_sd = 0
          )
        )
        fit <- nloptr(
          x0 = init, eval_f = objectivef,
          opts = control_restart
        )
        p()
        fit
      } else {
        list(status = fit$return_code, objective = Inf)
      }
    },
    future.seed = TRUE)
    
    logliks <- -unlist(lapply(out, "[[", "objective")) * n_obs
    return_codes <- unlist(lapply(out, "[[", "status"))
    successful <- which(return_codes > 0)
    if (length(successful) == 0) {
      warning_(
        c("All optimizations terminated due to error.",
          "Error of first restart: ", error_msg(return_codes[1]))
      )
    }
    optimum <- successful[which.max(logliks[successful])]
    pars <- out[[optimum]]$solution
    init <- list()
    init$pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
    init$A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
    if (C == 1L) {
      init$B <- create_eta_B_nhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_B
      )
    } else {
      init$B <- create_eta_multichannel_B_nhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_B
      )
    }
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
  }
  if (C == 1) {
    out <- EM_LBFGS_nhmm_singlechannel(
      init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
      Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
      n_obs, control$maxeval_em_dnm,
      control$ftol_abs, control$ftol_rel,
      control$xtol_abs, control$xtol_rel, 
      control$print_level, control_mstep$maxeval,
      control_mstep$ftol_abs, control_mstep$ftol_rel,
      control_mstep$xtol_abs, control_mstep$xtol_rel, 
      control_mstep$print_level, lambda, pseudocount)
  } else {
    out <- EM_LBFGS_nhmm_multichannel(
      init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
      Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
      n_obs, control$maxeval_em_dnm,
      control$ftol_abs, control$ftol_rel,
      control$xtol_abs, control$xtol_rel, 
      control$print_level, control_mstep$maxeval,
      control_mstep$ftol_abs, control_mstep$ftol_rel,
      control_mstep$xtol_abs, control_mstep$xtol_rel, 
      control_mstep$print_level, lambda, pseudocount)
  }
  if (out$return_code == 0) {
    init <- unlist(
      create_initial_values(
        stats::setNames(
          out[c("eta_pi", "eta_A", "eta_B")], c("pi", "A", "B")
        ), 
        model, 
        init_sd = 0
      )
    )
  } else {
    warning_(
      paste("EM-step terminated due to error:", error_msg(out$return_code),
            "Running DNM using initial values for EM.")
    )
    init <- unlist(create_initial_values(inits, model, init_sd))
  }
  out <- nloptr(
    x0 = init, eval_f = objectivef,
    opts = control
  )
  end_time <- proc.time()
  if (out$status < 0) {
    warning_(
      paste("Optimization terminated due to error:", error_msg(out$status))
    )
  }
  
  pars <- out$solution
  model$etas$pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$etas$A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  if (C == 1L) {
    model$etas$B <- create_eta_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B
    )
    model$gammas$B <- eta_to_gamma_cube(model$etas$B)
  } else {
    model$etas$B <- create_eta_multichannel_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B
    )
    model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
  }
  model$estimation_results <- list(
    loglik = -out$objective * n_obs, 
    return_code = out$status,
    message = out$message,
    iterations = out$iterations,
    logliks_of_restarts = if(restarts > 0L) logliks else NULL, 
    return_codes_of_restarts = if(restarts > 0L) return_codes else NULL,
    all_solutions = all_solutions,
    time = end_time - start_time,
    lambda = lambda,
    pseudocount = 0,
    method = "EM-DNM",
    algorithm = control$algorithm
  )
  model
}
