#' Estimate a Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_nhmm <- function(model, inits, init_sd, restarts, lambda, method, 
                     save_all_solutions = FALSE,
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
  
  if (isTRUE(control$maxeval < 0)) {
    pars <- unlist(create_initial_values(
      inits, S, M, init_sd, K_pi, K_A, K_B
    ))
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
    return(model)
  }
  all_solutions <- NULL
  if (method == "DNM") {
    need_grad <- grepl("NLOPT_LD_", control$algorithm)
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
      
      out <- future.apply::future_lapply(seq_len(restarts), function(i) {
        init <- unlist(create_initial_values(
          inits, S, M, init_sd, K_pi, K_A, K_B
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
      optimum <- successful[which.max(logliks[successful])]
      init <- out[[optimum]]$solution
      if (save_all_solutions) {
        all_solutions <- out
      }
    } else {
      init <- unlist(create_initial_values(
        inits, S, M, init_sd, K_pi, K_A, K_B
      ))
    }
    
    out <- nloptr(
      x0 = init, eval_f = objectivef,
      opts = control
    )
    end_time <- proc.time()
    if (out$status < 0) {
      warning_(paste("Optimization terminated due to error:", out$message))
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
      time = end_time - start_time
    )
  } else {
    start_time <- proc.time()
    if (restarts > 0L) {
      out <- future.apply::future_lapply(seq_len(restarts), function(i) {
        init <- create_initial_values(
          inits, S, M, init_sd, K_pi, K_A, K_B
        )
        if (C == 1) {
          EM_LBFGS_nhmm_singlechannel(
            init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
            Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
            n_obs, control_restart$maxeval,
            control_restart$ftol_abs, control_restart$ftol_rel,
            control_restart$xtol_abs, control_restart$xtol_rel, 
            control_restart$print_level, control_mstep$maxeval,
            control_mstep$ftol_abs, control_mstep$ftol_rel,
            control_mstep$xtol_abs, control_mstep$xtol_rel, 
            control_mstep$print_level, lambda)
        } else {
          EM_LBFGS_nhmm_multichannel(
            init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
            Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
            n_obs, control_restart$maxeval,
            control_restart$ftol_abs, control_restart$ftol_rel,
            control_restart$xtol_abs, control_restart$xtol_rel, 
            control_restart$print_level, control_mstep$maxeval,
            control_mstep$ftol_abs, control_mstep$ftol_rel,
            control_mstep$xtol_abs, control_mstep$xtol_rel, 
            control_mstep$print_level, lambda)
        }
      },
      future.seed = TRUE)
      
      logliks <- unlist(lapply(out, "[[", "penalized_logLik")) * n_obs
      optimum <- out[[which.max(logliks)]]
      init <- stats::setNames(
        optimum[c("eta_pi", "eta_A", "eta_B")], c("pi", "A", "B")
      )
      if (save_all_solutions) {
        all_solutions <- out
      }
    } else {
      init <- create_initial_values(
        inits, S, M, init_sd, K_pi, K_A, K_B
      )
    }
    if (C == 1) {
      out <- EM_LBFGS_nhmm_singlechannel(
        init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
        Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
        n_obs, control$maxeval,
        control$ftol_abs, control$ftol_rel,
        control$xtol_abs, control$xtol_rel, 
        control$print_level, control_mstep$maxeval,
        control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, lambda)
    } else {
      out <- EM_LBFGS_nhmm_multichannel(
        init$pi, model$X_pi, init$A, model$X_A, init$B, model$X_B, obs,
        Ti, icpt_only_pi, icpt_only_A, icpt_only_B, iv_A, iv_B, tv_A, tv_B,
        n_obs, control$maxeval,
        control$ftol_abs, control$ftol_rel,
        control$xtol_abs, control$xtol_rel, 
        control$print_level, control_mstep$maxeval,
        control_mstep$ftol_abs, control_mstep$ftol_rel,
        control_mstep$xtol_abs, control_mstep$xtol_rel, 
        control_mstep$print_level, lambda)
    }
    end_time <- proc.time()
    # if (out$status < 0) {
    #   warning_(paste("Optimization terminated due to error:", out$message))
    # }
    
    model$etas$pi[] <- out$eta_pi
    model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
    model$etas$A <- out$eta_A
    model$gammas$A <- eta_to_gamma_cube(model$etas$A)
    if (C == 1L) {
      model$etas$B[] <- out$eta_B
      model$gammas$B <- eta_to_gamma_cube(model$etas$B)
    } else {
      for (i in seq_len(C)) model$etas_B[[i]] <- out$eta_B[[i]]
      model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
    }
    
    model$estimation_results <- list(
      loglik = out$penalized_logLik * n_obs,
      iterations = out$iterations,
      logliks_of_restarts = if(restarts > 0L) logliks else NULL, 
      all_solutions = all_solutions,
      time = end_time - start_time,
      f_rel_change = out$relative_f_change,
      f_abs_change = out$absolute_f_change,
      x_rel_change = out$relative_x_change,
      x_abs_change = out$absolute_x_change
    ) 
  }
  model$estimation_results$lambda <- lambda
  model
}
