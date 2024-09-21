#' Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_mnhmm <- function(model, inits, init_sd, restarts, threads, hessian, ...) {
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  stopifnot_(
    checkmate::test_int(x = restarts, lower = 0L), 
    "Argument {.arg restarts} must be a single integer."
  )
  obs <- create_obsArray(model)
  if (model$n_channels == 1) {
    obs <- array(obs, dim(obs)[2:3])
  }
  M <- model$n_symbols
  S <- model$n_states
  D <- model$n_clusters
  T_ <- model$length_of_sequences
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
  
  n_i <- length(unlist(model$coefficients$gamma_pi_raw))
  n_s <- length(unlist(model$coefficients$gamma_A_raw))
  n_o <- length(unlist(model$coefficients$gamma_B_raw))
  n_d <- length(model$coefficients$gamma_omega_raw)
  X_i <- model$X_initial
  X_s <- model$X_transition
  X_o <- model$X_emission
  X_d <- model$X_cluster
  K_i <- nrow(X_i)
  K_s <- nrow(X_s)
  K_o <- nrow(X_o) 
  K_d <- nrow(X_d)
  
  dots <- list(...)
  if (is.null(dots$algorithm)) dots$algorithm <- "NLOPT_LD_LBFGS"
  need_grad <- grepl("NLOPT_LD_", dots$algorithm)
  
  if (model$n_channels == 1L) {
    if (need_grad) {
      objectivef <- function(pars) {
        gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        gamma_A_raw <- create_gamma_A_raw_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
        gamma_B_raw <- create_gamma_B_raw_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
        )
        gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- log_objective_mnhmm_singlechannel(
          gamma_pi_raw, X_i,
          gamma_A_raw, X_s,
          gamma_B_raw, X_o,
          gamma_omega_raw, X_d,
          obs)
        list(objective = - out$loglik,
             gradient = - unlist(out[-1]))
      }
    } else {
      objectivef <- function(pars) {
        gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        gamma_A_raw <- create_gamma_A_raw_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
        gamma_B_raw <- create_gamma_B_raw_mnhmm(
          pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
        )
        gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- forward_mnhmm_singlechannel(
          gamma_pi_raw, X_i,
          gamma_A_raw, X_s,
          gamma_B_raw, X_o,
          gamma_omega_raw, X_d,
          obs)
        
        - sum(apply(out[, T_, ], 2, logSumExp))
      }
    }
  } else {
    if (need_grad) {
      objectivef <- function(pars) {
        gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        gamma_A_raw <- create_gamma_A_raw_mnhmm(
          pars[n_i + seq_len(n_s)], 
          S, K_s, D
        )
        gamma_B_raw <- unlist(
          create_gamma_multichannel_B_raw_mnhmm(
            pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
          ), 
          recursive = FALSE
        )
        gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- log_objective_mnhmm_multichannel(
          gamma_pi_raw, X_i,
          gamma_A_raw, X_s,
          gamma_B_raw, X_o,
          gamma_omega_raw, X_d,
          obs, M)
        list(objective = - out$loglik,
             gradient = - unlist(out[-1]))
      }
    } else {
      objectivef <- function(pars) {
        gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
        gamma_A_raw <- create_gamma_A_raw_mnhmm(
          pars[n_i + seq_len(n_s)], S, K_s, D
          )
        gamma_B_raw <- unlist(
          create_gamma_multichannel_B_raw_mnhmm(
            pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
          ), 
          recursive = FALSE
        )
        gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
          pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
        )
        out <- forward_mnhmm_multichannel(
          gamma_pi_raw, X_i,
          gamma_A_raw, X_s,
          gamma_B_raw, X_o,
          gamma_omega_raw, X_d,
          obs, M)
        
        - sum(apply(out[, T_, ], 2, logSumExp))
      }
    }
  }
  
  if (restarts > 0L) {
    if (threads > 1L) {
      future::plan(future::multisession, workers = threads)
    } else {
      future::plan(future::sequential)
    }
    if (is.null(dots$maxeval)) dots$maxeval <- 1000L
    if (is.null(dots$print_level )) dots$print_level <- 0
    if (is.null(dots$xtol_rel)) dots$xtol_rel <- 1e-2
    if (is.null(dots$xtol_rel)) dots$ftol_rel <- 1e-4
    if (is.null(dots$check_derivatives)) dots$check_derivatives <- FALSE
    out <- future.apply::future_lapply(seq_len(restarts), function(i) {
      init <- unlist(create_initial_values(
        inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
      ))
      nloptr(
        x0 = init, eval_f = objectivef,
        opts = dots
      )
    },
    future.seed = TRUE)
    
    logliks <- unlist(lapply(out, "[[", "objective"))
    return_codes <- unlist(lapply(out, "[[", "status"))
    successful <- which(return_codes > 0)
    optimum <- successful[which.max(logliks[successful])]
    init <- out[[optimum]]$solution
  } else {
    init <- unlist(create_initial_values(
      inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
    ))
  }
  dots <- list(...)
  if (is.null(dots$algorithm)) dots$algorithm <- "NLOPT_LD_LBFGS"
  if (is.null(dots$maxeval)) dots$maxeval <- 10000L
  if (is.null(dots$xtol_rel)) dots$xtol_rel <- 1e-6
  if (is.null(dots$xtol_rel)) dots$ftol_rel <- 1e-12
  if (is.null(dots$check_derivatives)) dots$check_derivatives <- FALSE
  out <- nloptr(
    x0 = init, eval_f = objectivef,
    opts = dots
  )
  if (out$status < 0) {
    warning_(paste("Optimization terminated due to error:", out$message))
  }
  pars <- out$solution
  model$coefficients$gamma_pi_raw <- create_gamma_pi_raw_mnhmm(
    pars[seq_len(n_i)], S, K_i, D
    )
  model$coefficients$gamma_A_raw <- create_gamma_A_raw_mnhmm(
    pars[n_i + seq_len(n_s)], S, K_s, D
    )
  if (model$n_channels == 1L) {
    model$coefficients$gamma_B_raw <- create_gamma_B_raw_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
  } else {
    model$coefficients$gamma_B_raw <- create_gamma_multichannel_B_raw_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
  }
  model$coefficients$gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
    pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
  )
  
  if (!isFALSE(hessian)) {
    if (isTRUE(hessian)){
      hessian <- numDeriv::jacobian(objectivef, pars)
    } else {
      hessian <- do.call(
        numDeriv::jacobian,
        c(hessian, list(func = objectivef, x = pars))
      )
    }
  } else {
    hessian <- NULL
  }
  
  model$estimation_results <- list(
    hessian = hessian,
    loglik = out$objective, 
    return_code = out$status,
    message = out$message,
    iterations = out$iterations,
    logliks_of_restarts = if(restarts > 0L) logliks else NULL, 
    return_codes_of_restarts = if(restarts > 0L) return_codes else NULL
  )
  model
}
