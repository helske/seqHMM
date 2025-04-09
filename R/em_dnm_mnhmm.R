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
  K_omega <- nrow(X_omega)
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- vapply(X_B, \(x) nrow(x), 1L)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  need_grad <- grepl("NLOPT_LD_", control$algorithm)
  obs <- create_obsArray(model)
  all_solutions <- NULL
  if (need_grad) {
    objectivef <- function(pars) {
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
      out <- log_objective_mnhmm(
        obs, Ti, M, X_pi, X_A, X_B, X_omega,
        icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
        iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega
      )
      list(
        objective = - (out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs, 
        gradient = - (unlist(out[-1]) - lambda * pars) / n_obs
      )
    }
  } else {
    objectivef <- function(pars) {
      eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
      eta_A <- create_eta_A_mnhmm(
        pars[np_pi + seq_len(np_A)], S, K_A, D
      )
      eta_B <- create_eta_B_mnhmm(
        pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
      )
      eta_omega <- create_eta_omega_mnhmm(
        pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
      )
      out <- forward_mnhmm(
        obs, Ti, M, X_pi, X_A, X_B, X_omega,
        icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
        iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega
      )
      
      - (sum(apply(out[, T_, ], 2, logSumExp)) - 0.5 * lambda * sum(pars^2)) / n_obs
    }
  }
  
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
      em_return_code <- fit$return_code
      if (em_return_code >= 0) {
        init <- unlist(fit[c("eta_omega", "eta_pi", "eta_A", "eta_B")])
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
        p()
        list(status = -1000 + em_return_code, objective = NaN)
      }
    }
    p <- progressr::progressor(along = seq_len(restarts))
    original_options <- options(future.globals.maxSize = Inf)
    on.exit(options(original_options))
    base_init <- create_initial_values(inits, model, init_sd = 0)
    u <- t(
      stats::qnorm(
        lhs::maximinLHS(restarts, length(unlist(base_init))), 
        sd = init_sd
      )
    )
    out <- future.apply::future_lapply(
      seq_len(restarts), \(i) .fun(i, base_init, u[, i]), future.seed = TRUE
    )
    logliks <- -unlist(lapply(out, "[[", "objective")) * n_obs
    return_codes <- unlist(lapply(out, "[[", "status"))
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
      pars <- out[[optimum]]$solution
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
      all_solutions <- out
    }
  } else {
    init <- create_initial_values(inits, model, init_sd)
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
    em_return_code <- out$return_code
    if (em_return_code >= 0) {
      init <- out[c("eta_omega", "eta_pi", "eta_A", "eta_B")]
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
    loglik <- NaN
  } else {
    loglik <- -out$objective * n_obs
  }
  
  pars <- out$solution
  model$etas$pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
  model$gammas$pi <- drop(eta_to_gamma_mat_field(model$etas$pi))
  model$etas$A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
  model$gammas$A <- drop(eta_to_gamma_cube_field(model$etas$A))
  model$etas$B <- create_eta_B_mnhmm(
    pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
  )
  model$gammas$B <- split(
    eta_to_gamma_cube_2d_field(model$etas$B), seq_len(D)
  )
  model$etas$omega <- create_eta_omega_mnhmm(
    pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
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
