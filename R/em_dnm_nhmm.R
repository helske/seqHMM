em_dnm_nhmm <- function(model, inits, init_sd, restarts, lambda,
                        bound, control, control_restart, control_mstep, 
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
      if (C == 1) {
        fit <- EM_LBFGS_nhmm_singlechannel(
          init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, obs,
          Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
          n_obs, control$maxeval_em_dnm,
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
          n_obs, control$maxeval_em_dnm,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, lambda, bound)
      }
      em_return_code <- fit$return_code
      if (em_return_code >= 0) {
        init <- unlist(
          create_initial_values(
            fit[c("eta_pi", "eta_A", "eta_B")], 
            model, 
            init_sd = 0
          )
        )
        fit <- nloptr(
          x0 = init, eval_f = objectivef, lb = -rep(bound, length(init)), 
          ub = rep(bound, length(init)), opts = control_restart
        )
        if (fit$status == -1) {
          grad_norm <- if (need_grad) {
            max(abs(objectivef(fit$solution)$gradient))
          } else {
            1
          }
          relative_change <- abs(fit$objective - objectivef(init)$objective) / 
            (abs(fit$objective) + 1e-12)
          if ((grad_norm < 1e-8 && is.finite(fit$objective)) || 
              relative_change < control_restart$ftol_rel) {
            fit$status <- 7
          }
        }
        p()
        fit
      } else {
        list(status = -1000 + em_return_code, objective = NaN)
      }
    },
    future.seed = TRUE)
    
    logliks <- -unlist(lapply(out, "[[", "objective")) * n_obs
    return_codes <- unlist(lapply(out, "[[", "status"))
    successful <- which(return_codes > 0)
    if (length(successful) == 0) {
      warning_(
        c("All restarts terminated due to error.",
          "Error of first restart: ", return_msg(return_codes[1]),
          "Running DNM using initial values for EM.")
      )
      em_return_code <- return_codes[1] + 1000
      init <- create_initial_values(inits, model, init_sd)
    } else {
      em_return_code <- 0 # generic success
      optimum <- successful[which.max(logliks[successful])]
      pars <- out[[optimum]]$solution
      init <- list()
      init$eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
      init$eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
      if (C == 1L) {
        init$eta_B <- create_eta_B_nhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B
        )
      } else {
        init$eta_B <- create_eta_multichannel_B_nhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B
        )
      }
    }
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
    if (C == 1) {
      out <- EM_LBFGS_nhmm_singlechannel(
        init$eta_pi, model$X_pi, init$eta_A, model$X_A, init$eta_B, model$X_B, obs,
        Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
        n_obs, control$maxeval_em_dnm,
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
        n_obs, control$maxeval_em_dnm,
        control$ftol_abs, control$ftol_rel,
        control$xtol_abs, control$xtol_rel, 
        control$print_level, control_mstep$maxeval,
        control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, lambda, bound)
    }
    em_return_code <- out$return_code
    if (em_return_code >= 0) {
      init <- out[c("eta_pi", "eta_A", "eta_B")]
    } else {
      warning_(
        paste("EM-step terminated due to error:", return_msg(em_return_code),
              "Running DNM using initial values for EM.")
      )
      init <- create_initial_values(inits, model, init_sd)
    }
  }
  out <- nloptr(
    x0 = unlist(init), eval_f = objectivef, 
    lb = -rep(bound, length(unlist(init))), 
    ub = rep(bound, length(unlist(init))), opts = control
  )
  if (out$status == -1) {
    grad_norm <- if (need_grad) {
      max(abs(objectivef(out$solution)$gradient))
    } else {
      1
    }
    relative_change <- abs(out$objective - objectivef(unlist(init))$objective) / 
      (abs(out$objective) + 1e-12)
    if ((grad_norm < 1e-8 && is.finite(out$objective)) || 
        relative_change < control$ftol_rel) {
      out$status <- 7
    }
  }
  if (out$status < 0) {
    warning_(
      paste("Optimization terminated due to error:", return_msg(out$status))
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
    lambda = lambda,
    bound = bound,
    method = "EM-DNM",
    algorithm = control$algorithm,
    EM_return_code = em_return_code
  )
  model
}
