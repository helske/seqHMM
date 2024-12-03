em_dnm_mnhmm <- function(model, inits, init_sd, restarts, lambda, 
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
  need_grad <- grepl("NLOPT_LD_", control$algorithm)
  obs <- create_obsArray(model)
  if (C == 1) {
    obs <- array(obs, dim(obs)[2:3])
  }
  all_solutions <- NULL
  if (C == 1L) {
    if (need_grad) {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
        eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
        eta_B <- create_eta_B_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
        )
        out <- log_objective_mnhmm_singlechannel(
          eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B,
          obs, Ti, icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B, 
          iv_A, iv_B, tv_A, tv_B
        )
        list(
          objective = - (out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs, 
          gradient = - (unlist(out[-1]) - lambda * pars) / n_obs
        )
      }
    } else {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
        eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
        eta_B <- create_eta_B_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
        )
        out <- forward_mnhmm_singlechannel(
          eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B, 
          obs, Ti, icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B, 
          iv_A, iv_B, tv_A, tv_B
        )
        
        - (sum(apply(out[, T_, ], 2, logSumExp)) - 0.5 * lambda * sum(pars^2)) / n_obs
      }
    }
  } else {
    if (need_grad) {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
        eta_A <- create_eta_A_mnhmm(
          pars[n_i + seq_len(n_s)], 
          S, K_A, D
        )
        eta_B <- unlist(
          create_eta_multichannel_B_mnhmm(
            pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
          ),
          recursive = FALSE
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
        )
        out <- log_objective_mnhmm_multichannel(
          eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B,
          obs, Ti, icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B, 
          iv_A, iv_B, tv_A, tv_B
        )
        list(
          objective = - (out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs, 
          gradient = - (unlist(out[-1]) - lambda * pars) / n_obs
        )
      }
    } else {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
        eta_A <- create_eta_A_mnhmm(
          pars[n_i + seq_len(n_s)], S, K_A, D
        )
        eta_B <- unlist(
          create_eta_multichannel_B_mnhmm(
            pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
          ),
          recursive = FALSE
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
        )
        out <- forward_mnhmm_multichannel(
          eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B, 
          obs, Ti, icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B, 
          iv_A, iv_B, tv_A, tv_B
        )
        
        - (sum(apply(out[, T_, ], 2, logSumExp)) - 0.5 * lambda * sum(pars^2)) / n_obs
      }
    }
  }
  if (restarts > 0L) {
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    out <- future.apply::future_lapply(seq_len(restarts), function(i) {
      init <- create_initial_values(inits, model, init_sd)
      if (C == 1) {
        fit <- EM_LBFGS_mnhmm_singlechannel(
          init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
          init$B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
          icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
          n_obs, control_restart$maxeval,
          control_restart$ftol_abs, control_restart$ftol_rel,
          control_restart$xtol_abs, control_restart$xtol_rel, 
          control_restart$print_level, control_mstep$maxeval,
          control_mstep$ftol_abs, control_mstep$ftol_rel,
          control_mstep$xtol_abs, control_mstep$xtol_rel, 
          control_mstep$print_level, lambda, bound)
      } else {
        eta_B <- unlist(init$B, recursive = FALSE)
        fit <- EM_LBFGS_mnhmm_multichannel(
          init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
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
      em_return_code <- fit$return_code
      if (em_return_code >= 0) {
        init <- unlist(
          create_initial_values(
            stats::setNames(
              fit[c("eta_omega", "eta_pi", "eta_A", "eta_B")], 
              c("omega", "pi", "A", "B")
            ), 
            model, 
            init_sd = 0
          )
        )
        fit <- nloptr(
          x0 = init, eval_f = objectivef, lb = -rep(bound, length(init)), 
          ub = rep(bound, length(init)), opts = control_restart
        )
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
        c("All optimizations terminated due to error.",
          "Error of first restart: ", error_msg(return_codes[1]),
          "Running DNM using initial values for EM.")
      )
      init <- create_initial_values(inits, model, init_sd)
      em_return_code <- return_codes[1] + 1000
    } else {
      em_return_code <- 0 # generic success
      optimum <- successful[which.max(logliks[successful])]
      pars <- out[[optimum]]$solution
      init <- list()
      init$pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
      init$A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
      if (C == 1L) {
        init$B <- create_eta_B_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
        )
      } else {
        init$B <- create_eta_multichannel_B_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
        )
      }
      init$omega <- create_eta_omega_mnhmm(
        pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
      )
    }
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
    if (C == 1) {
      out <- EM_LBFGS_mnhmm_singlechannel(
        init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
        init$B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
        icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
        n_obs, control$maxeval,
        control$ftol_abs, control$ftol_rel,
        control$xtol_abs, control$xtol_rel, 
        control$print_level, control_mstep$maxeval,
        control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, lambda, bound)
    } else {
      eta_B <- unlist(init$B, recursive = FALSE)
      out <- EM_LBFGS_mnhmm_multichannel(
        init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
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
    em_return_code <- out$return_code
    if (em_return_code >= 0) {
      init <- create_initial_values(
        stats::setNames(
          out[c("eta_omega", "eta_pi", "eta_A", "eta_B")], 
          c("omega", "pi", "A", "B")
        ), 
        model, 
        init_sd = 0
      )
    } else {
      warning_(
        paste("EM-step terminated due to error:", error_msg(em_return_code),
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
  if (out$status < 0) {
    warning_(
      paste("Optimization terminated due to error:", error_msg(out$status))
    )
    loglik <- NaN
  } else {
    loglik <- -out$objective * n_obs
  }
  
  pars <- out$solution
  model$etas$pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
  model$gammas$pi <- c(eta_to_gamma_mat_field(model$etas$pi))
  model$etas$A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
  model$gammas$A <- c(eta_to_gamma_cube_field(model$etas$A))
  if (C == 1L) {
    model$etas$B <- create_eta_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
    )
    model$gammas$B <- c(eta_to_gamma_cube_field(model$etas$B))
  } else {
    model$etas$B <- create_eta_multichannel_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
    )
    l <- lengths(model$etas$B)
    gamma_B <- c(eta_to_gamma_cube_field(unlist(model$etas$B, recursive = FALSE)))
    model$gammas$B <- unname(split(gamma_B, rep(seq_along(l), l)))
  }
  model$etas$omega <- create_eta_omega_mnhmm(
    pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
  )
  model$gammas$omega <- eta_to_gamma_mat(model$etas$omega)
  model$estimation_results <- list(
    loglik = loglik, 
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
