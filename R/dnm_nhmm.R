dnm_nhmm <- function(model, inits, init_sd, restarts, lambda, bound, control, 
                     control_restart, save_all_solutions) {
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
      eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
      eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
      eta_B <- create_eta_B_nhmm(
        pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
      )
      out <- log_objective_nhmm(
        obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B,
        iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B
      )
      list(
        objective = - (out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs, 
        gradient = - (unlist(out[-1]) - lambda * pars) / n_obs
      )
    }
  } else {
    objectivef <- function(pars) {
      eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
      eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
      eta_B <- create_eta_B_nhmm(
        pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
      )
      out <- forward_nhmm(
        obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi,icpt_only_A, icpt_only_B,
        iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B
      )
      - (sum(apply(out[, T_, ], 2, logSumExp)) - 0.5 * lambda * sum(pars^2)) / n_obs
    }
  }
  if (restarts > 0L) {
    .fun <- function(i, u, base_init) {
      init <- base_init + u
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
          "Trying once more. ")
      )
      init <- unlist(create_initial_values(inits, model, init_sd))
    } else {
      optimum <- successful[which.max(logliks[successful])]
      init <- out[[optimum]]$solution
    }
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- unlist(create_initial_values(inits, model, init_sd))
  }
  
  out <- nloptr(
    x0 = init, eval_f = objectivef, lb = -rep(bound, length(init)), 
    ub = rep(bound, length(init)), opts = control
  )
  if (out$status == -1) {
    grad_norm <- if (need_grad) {
      max(abs(objectivef(out$solution)$gradient))
    } else {
      1
    }
    relative_change <- abs(out$objective - objectivef(init)$objective) / 
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
  model$etas$pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
  model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
  model$etas$A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
  model$gammas$A <- eta_to_gamma_cube(model$etas$A)
  model$etas$B <- create_eta_B_nhmm(
    pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
  )
  model$gammas$B <- drop(eta_to_gamma_cube_field(model$etas$B))
  
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
    method = "DNM", 
    algorithm = control$algorithm
  )
  model
}
