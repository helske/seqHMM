#' Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_mnhmm <- function(model, inits, init_sd, restarts, lambda, method,
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
  iv_pi <- attr(model$X_pi, "iv")
  iv_A <- attr(model$X_A, "iv")
  iv_B <- attr(model$X_B, "iv")
  iv_omega <- attr(model$X_omega, "iv")
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
      model$gammas$B <- split(gamma_B, rep(seq_along(l), l))
    }
    model$etas$omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
    )
    model$gammas$omega <- eta_to_gamma_mat(
      model$etas$omega
    )
    return(model)
  }
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
          obs, iv_omega, iv_pi, iv_A, iv_B, tv_A, tv_B, Ti
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
          obs, iv_omega, iv_pi, iv_A, iv_B, tv_A, tv_B, Ti
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
          obs, iv_omega, iv_pi, iv_A, iv_B, tv_A, tv_B, Ti
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
          eta_pi, X_pi,
          eta_A, X_A,
          eta_B, X_B,
          eta_omega, X_omega,
          obs, M)

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
    warning_(paste("Optimization terminated due to error:", out$message))
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
    model$gammas$B <- split(gamma_B, rep(seq_along(l), l))
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
    time = end_time - start_time
  )
  model
}
