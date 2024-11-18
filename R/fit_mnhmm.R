#' Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_mnhmm <- function(model, inits, init_sd, restarts, lambda, method, 
                      pseudocount, save_all_solutions = FALSE, 
                      control_restart = list(), control_mstep = list(), ...) {
  
  stopifnot_(
    checkmate::test_int(x = restarts, lower = 0L), 
    "Argument {.arg restarts} must be a single integer."
  )
  control <- utils::modifyList(
    list(
      ftol_abs = 1e-8,
      ftol_rel = 1e-8,
      xtol_abs = 1e-4,
      xtol_rel = 1e-4,
      maxeval = 1e4,
      print_level = 0,
      algorithm = "NLOPT_LD_LBFGS"
    ),
    list(...)
  )
  control_restart <- utils::modifyList(control, control_restart)
  control_mstep <- utils::modifyList(control, control_mstep)
  
  M <- model$n_symbols
  S <- model$n_states
  D <- model$n_clusters
  T_ <- model$length_of_sequences
  C <- model$n_channels
  obs <- create_obsArray(model)
  if (C == 1) {
    obs <- array(obs, dim(obs)[2:3])
  }
  if (identical(inits, "random")) {
    inits <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL)
  } else {
    if (is.null(inits$initial_probs)) inits$initial_probs <- NULL
    if (is.null(inits$transition_probs)) inits$transition_probs <- NULL
    if (is.null(inits$emission_probs)) inits$emission_probs <- NULL
  }
  
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
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  X_omega <- model$X_omega
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B) 
  K_omega <- nrow(X_omega)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  
  if (isTRUE(control$maxeval < 0)) {
    pars <- unlist(create_initial_values(
      inits, S, M, init_sd, K_pi, K_A, K_B, K_omega, D
    ))
    model$etas$pi <- create_eta_pi_mnhmm(
      pars[seq_len(n_i)], S, K_pi, D
    )
    model$gammas$pi <- c(eta_to_gamma_mat_field(
      model$etas$pi
    ))
    model$etas$A <- create_eta_A_mnhmm(
      pars[n_i + seq_len(n_s)], S, K_A, D
    )
    model$gammas$A <- c(eta_to_gamma_cube_field(
      model$etas$A
    ))
    if (C == 1L) {
      model$etas$B <- create_eta_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
      )
      model$gammas$B <- c(eta_to_gamma_cube_field(
        model$etas$B
      ))
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
    model$gammas$omega <- eta_to_gamma_mat(
      model$etas$omega
    )
    return(model)
  }
  all_solutions <- NULL
  if (method == "DNM") {
    need_grad <- grepl("NLOPT_LD_", control$algorithm)
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
    all_solutions <- NULL
    start_time <- proc.time()
    if (restarts > 0L) {
      out <- future.apply::future_lapply(seq_len(restarts), function(i) {
        init <- unlist(create_initial_values(
          inits, S, M, init_sd, K_pi, K_A, K_B, K_omega, D
        ))
        nloptr(
          x0 = init, eval_f = objectivef,
          opts = control_restart
        )
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
      init <- out[[optimum]]$solution
      if (save_all_solutions) {
        all_solutions <- out
      }
    } else {
      init <- unlist(create_initial_values(
        inits, S, M, init_sd, K_pi, K_A, K_B, K_omega, D
      ))
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
    model$etas$pi <- create_eta_pi_mnhmm(
      pars[seq_len(n_i)], S, K_pi, D
    )
    model$gammas$pi <- c(eta_to_gamma_mat_field(
      model$etas$pi
    ))
    model$etas$A <- create_eta_A_mnhmm(
      pars[n_i + seq_len(n_s)], S, K_A, D
    )
    model$gammas$A <- c(eta_to_gamma_cube_field(
      model$etas$A
    ))
    if (C == 1L) {
      model$etas$B <- create_eta_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
      )
      model$gammas$B <- c(eta_to_gamma_cube_field(
        model$etas$B
      ))
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
    model$gammas$omega <- eta_to_gamma_mat(
      model$etas$omega
    )
    
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
      pseudocount = pseudocount,
      method = method
    )
    return(model)
  }
  if (method == "EM") {
    start_time <- proc.time()
    if (restarts > 0L) {
      out <- future.apply::future_lapply(seq_len(restarts), function(i) {
        init <- create_initial_values(
          inits, S, M, init_sd, K_pi, K_A, K_B, K_omega, D
        )
        if (C == 1) {
          EM_LBFGS_mnhmm_singlechannel(
            init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
            init$B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
            icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
            n_obs, control_restart$maxeval,
            control_restart$ftol_abs, control_restart$ftol_rel,
            control_restart$xtol_abs, control_restart$xtol_rel, 
            control_restart$print_level, control_mstep$maxeval,
            control_mstep$ftol_abs, control_mstep$ftol_rel,
            control_mstep$xtol_abs, control_mstep$xtol_rel, 
            control_mstep$print_level, lambda, pseudocount)
        } else {
          EM_LBFGS_mnhmm_multichannel(
            init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
            init$B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
            icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
            n_obs, control_restart$maxeval,
            control_restart$ftol_abs, control_restart$ftol_rel,
            control_restart$xtol_abs, control_restart$xtol_rel, 
            control_restart$print_level, control_mstep$maxeval,
            control_mstep$ftol_abs, control_mstep$ftol_rel,
            control_mstep$xtol_abs, control_mstep$xtol_rel, 
            control_mstep$print_level, lambda, pseudocount)
        }
      },
      future.seed = TRUE)
      return_codes <- unlist(lapply(out, "[[", "return_code"))
      if (all(return_codes < 0)) {
        warning_(
          c("All optimizations terminated due to error.",
            "Error of first restart: ", error_msg(return_codes[1]))
        )
      }
      logliks <- unlist(lapply(out, "[[", "penalized_logLik")) * n_obs
      optimum <- out[[which.max(logliks)]]
      init <- stats::setNames(
        optimum[c("eta_omega", "eta_pi", "eta_A", "eta_B")], 
        c("omega", "pi", "A", "B")
      )
      if (save_all_solutions) {
        all_solutions <- out
      }
    } else {
      init <- create_initial_values(
        inits, S, M, init_sd, K_pi, K_A, K_B, K_omega, D
      )
    }
    if (C == 1) {
      out <- EM_LBFGS_mnhmm_singlechannel(
        init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
        init$B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
        icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
        n_obs, control_restart$maxeval,
        control_restart$ftol_abs, control_restart$ftol_rel,
        control_restart$xtol_abs, control_restart$xtol_rel, 
        control_restart$print_level, control_mstep$maxeval,
        control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, lambda, pseudocount)
    } else {
      eta_B <- unlist(init$B, recursive = FALSE)
      out <- EM_LBFGS_mnhmm_multichannel(
        init$omega, model$X_omega, init$pi, model$X_pi, init$A, model$X_A, 
        eta_B, model$X_B, obs, Ti, icpt_only_omega, icpt_only_pi, 
        icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
        n_obs, control_restart$maxeval,
        control_restart$ftol_abs, control_restart$ftol_rel,
        control_restart$xtol_abs, control_restart$xtol_rel, 
        control_restart$print_level, control_mstep$maxeval,
        control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, lambda, pseudocount)
    }
    end_time <- proc.time()
    if (out$return_code < 0) {
      warning_(
        paste("Optimization terminated due to error:", error_msg(out$return_code))
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
      loglik = out$penalized_logLik,
      penalty = out$penalty_term,
      iterations = out$iterations,
      return_code = out$return_code,
      logliks_of_restarts = if(restarts > 0L) logliks else NULL, 
      return_codes_of_restarts = if(restarts > 0L) return_codes else NULL,
      all_solutions = all_solutions,
      time = end_time - start_time,
      f_rel_change = out$relative_f_change,
      f_abs_change = out$absolute_f_change,
      x_rel_change = out$relative_x_change,
      x_abs_change = out$absolute_x_change,
      lambda = lambda,
      pseudocount = pseudocount,
      method = method
    )
    return(model)
  }
}
