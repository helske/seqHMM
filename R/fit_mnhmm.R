#' Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_mnhmm <- function(model, inits, init_sd, restarts, threads, 
                      save_all_solutions = FALSE, ...) {
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  stopifnot_(
    checkmate::test_int(x = restarts, lower = 0L), 
    "Argument {.arg restarts} must be a single integer."
  )
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
  iv_pi <- attr(model$X_initial, "iv")
  iv_A <- attr(model$X_transition, "iv")
  iv_B <- attr(model$X_emission, "iv")
  iv_omega <- attr(model$X_cluster, "iv")
  tv_A <- attr(model$X_transition, "tv")
  tv_B <- attr(model$X_emission, "tv")
  X_i <- model$X_initial
  X_s <- model$X_transition
  X_o <- model$X_emission
  X_d <- model$X_cluster
  K_i <- nrow(X_i)
  K_s <- nrow(X_s)
  K_o <- nrow(X_o) 
  K_d <- nrow(X_d)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  dots <- list(...)
  
  if (isTRUE(dots$maxeval < 0)) {
    pars <- unlist(create_initial_values(
      inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
    ))
    model$etas$pi <- create_eta_pi_mnhmm(
      pars[seq_len(n_i)], S, K_i, D
    )
    model$gammas$pi <- c(eta_to_gamma_mat_field(
      model$etas$pi
    ))
    model$etas$A <- create_eta_A_mnhmm(
      pars[n_i + seq_len(n_s)], S, K_s, D
    )
    model$gammas$A <- c(eta_to_gamma_cube_field(
      model$etas$A
    ))
    if (C == 1L) {
      model$etas$B <- create_eta_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      )
      model$gammas$B <- c(eta_to_gamma_cube_field(
        model$etas$B
      ))
    } else {
      model$etas$B <- create_eta_multichannel_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      )
      l <- lengths(model$etas$B)
      gamma_B <- c(eta_to_gamma_cube_field(unlist(model$etas$B, recursive = FALSE)))
      model$gammas$B <- split(gamma_B, rep(seq_along(l), l))
    }
    model$etas$omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    model$gammas$omega <- eta_to_gamma_mat(
      model$etas$omega
    )
    return(model)
  }
  if (is.null(dots$algorithm)) 
    dots$algorithm <- "NLOPT_LD_LBFGS"
  need_grad <- grepl("NLOPT_LD_", dots$algorithm)
  if (is.null(dots$maxeval)) 
    dots$maxeval <- 10000L
  if (is.null(dots$xtol_abs)) 
    dots$xtol_abs <- rep(1e-8, attr(model, "df"))
  if (is.null(dots$xtol_rel)) 
    dots$xtol_rel <- 0
  if (is.null(dots$ftol_abs)) 
    dots$ftol_abs <- 1e-8
  if (is.null(dots$ftol_rel)) 
    dots$ftol_rel <- 1e-8
  if (is.null(dots$check_derivatives)) 
    dots$check_derivatives <- FALSE
  
  if (C == 1L) {
    Qs <- t(create_Q(S))
    Qm <- t(create_Q(M))
    Qd <- t(create_Q(D))
    if (need_grad) {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
        eta_B <- create_eta_B_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- log_objective_mnhmm_singlechannel(
          Qs, Qm, Qd, eta_pi, X_i, eta_A, X_s, eta_B, X_o,
          eta_omega, X_d, obs, iv_pi, iv_A, iv_B, tv_A, tv_B, iv_omega, 
          Ti
        )
        list(
          objective = - out$loglik / n_obs, 
          gradient = - unlist(out[-1]) / n_obs
        )
      }
    } else {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
        eta_B <- create_eta_B_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- forward_mnhmm_singlechannel(
          eta_pi, X_i, eta_A, X_s, eta_B, X_o, 
          eta_omega, X_d, obs
        )
        
        - sum(apply(out[, T_, ], 2, logSumExp)) / n_obs
      }
    }
  } else {
    Qs <- t(create_Q(S))
    Qm <- lapply(M, function(m) t(create_Q(m)))
    Qd <- t(create_Q(D))
    if (need_grad) {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        eta_A <- create_eta_A_mnhmm(
          pars[n_i + seq_len(n_s)], 
          S, K_s, D
        )
        eta_B <- unlist(
          create_eta_multichannel_B_mnhmm(
            pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
          ),
          recursive = FALSE
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- log_objective_mnhmm_multichannel(
          Qs, Qm, Qd, eta_pi, X_i, eta_A, X_s, eta_B, X_o,
          eta_omega, X_d, obs, M, iv_pi, iv_A, iv_B, tv_A, tv_B, iv_omega,
          Ti
        )
        list(
          objective = - out$loglik / n_obs,
          gradient = - unlist(out[-1]) / n_obs
        )
      }
    } else {
      objectivef <- function(pars) {
        eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        eta_A <- create_eta_A_mnhmm(
          pars[n_i + seq_len(n_s)], S, K_s, D
        )
        eta_B <- unlist(
          create_eta_multichannel_B_mnhmm(
            pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
          ),
          recursive = FALSE
        )
        eta_omega <- create_eta_omega_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- forward_mnhmm_multichannel(
          eta_pi, X_i,
          eta_A, X_s,
          eta_B, X_o,
          eta_omega, X_d,
          obs, M)

        - sum(apply(out[, T_, ], 2, logSumExp)) / n_obs
      }
    }
  }
  all_solutions <- NULL
  start_time <- proc.time()
  if (restarts > 0L) {
    if (threads > 1L) {
      future::plan(future::multisession, workers = threads)
    } else {
      future::plan(future::sequential)
    }
    dots$control_restart$algorithm <- dots$algorithm
    if (is.null(dots$control_restart$maxeval)) 
      dots$control_restart$maxeval <- dots$maxeval
    if (is.null(dots$control_restart$print_level)) 
      dots$control_restart$print_level <- 0
    if (is.null(dots$control_restart$xtol_abs)) 
      dots$control_restart$xtol_abs <-dots$xtol_abs
    if (is.null(dots$control_restart$ftol_abs)) 
      dots$control_restart$ftol_abs <-  dots$ftol_abs
    if (is.null(dots$control_restart$xtol_rel)) 
      dots$control_restart$xtol_rel <- dots$xtol_rel
    if (is.null(dots$control_restart$ftol_rel)) 
      dots$control_restart$ftol_rel <- dots$ftol_rel
    out <- future.apply::future_lapply(seq_len(restarts), function(i) {
      init <- unlist(create_initial_values(
        inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
      ))
      nloptr(
        x0 = init, eval_f = objectivef,
        opts = dots$control_restart
      )
    },
    future.seed = TRUE)
    
    logliks <- -unlist(lapply(out, "[[", "objective"))
    return_codes <- unlist(lapply(out, "[[", "status"))
    successful <- which(return_codes > 0)
    optimum <- successful[which.max(logliks[successful])]
    init <- out[[optimum]]$solution
    if (save_all_solutions) {
      all_solutions <- out
    }
  } else {
    init <- unlist(create_initial_values(
      inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
    ))
  }
  out <- nloptr(
    x0 = init, eval_f = objectivef,
    opts = dots
  )
  end_time <- proc.time()
  if (out$status < 0) {
    warning_(paste("Optimization terminated due to error:", out$message))
  }
  pars <- out$solution
  model$etas$pi <- create_eta_pi_mnhmm(
    pars[seq_len(n_i)], S, K_i, D
  )
  model$gammas$pi <- c(eta_to_gamma_mat_field(
    model$etas$pi
  ))
  model$etas$A <- create_eta_A_mnhmm(
    pars[n_i + seq_len(n_s)], S, K_s, D
  )
  model$gammas$A <- c(eta_to_gamma_cube_field(
    model$etas$A
  ))
  if (C == 1L) {
    model$etas$B <- create_eta_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
    model$gammas$B <- c(eta_to_gamma_cube_field(
      model$etas$B
    ))
  } else {
    model$etas$B <- create_eta_multichannel_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
    l <- lengths(model$etas$B)
    gamma_B <- c(eta_to_gamma_cube_field(unlist(model$etas$B, recursive = FALSE)))
    model$gammas$B <- split(gamma_B, rep(seq_along(l), l))
  }
  model$etas$omega <- create_eta_omega_mnhmm(
    pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
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
