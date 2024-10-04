#' Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_mnhmm <- function(model, inits, init_sd, restarts, threads, penalty, ...) {
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
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  n_d <- attr(model, "np_omega")
  iv_pi <- attr(model, "iv_pi")
  iv_A <- attr(model, "iv_A")
  iv_B <- attr(model, "iv_B")
  iv_omega <- attr(model, "iv_omega")
  tv_A <- attr(model, "tv_A")
  tv_B <- attr(model, "tv_B")
  X_i <- model$X_initial
  X_s <- model$X_transition
  X_o <- model$X_emission
  X_d <- model$X_cluster
  K_i <- nrow(X_i)
  K_s <- nrow(X_s)
  K_o <- nrow(X_o) 
  K_d <- nrow(X_d)
  Ti <- model$sequence_lengths
  
  dots <- list(...)
  
  if (isTRUE(dots$maxeval < 0)) {
    pars <- unlist(create_initial_values(
      inits, S, M, init_sd, K_i, K_s, K_o, K_d, D
    ))
    model$etas$pi <- create_eta_pi_mnhmm(
      pars[seq_len(n_i)], S, K_i, D
    )
    model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
    model$etas$A <- create_eta_A_mnhmm(
      pars[n_i + seq_len(n_s)], S, K_s, D
    )
    model$gammas$A <- eta_to_gamma_cube(model$etas$A)
    if (model$n_channels == 1L) {
      model$etas$B <- create_eta_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      )
      model$gammas$B <- eta_to_gamma_cube(model$etas$B)
    } else {
      model$etas$B <- create_eta_multichannel_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      )

      model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
    }
    model$etas$omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    model$gammas$omega <- eta_to_gamma_mat(model$etas$omega)
    return(model)
  }
  
  if (is.null(dots$algorithm)) dots$algorithm <- "NLOPT_LD_LBFGS"
  need_grad <- grepl("NLOPT_LD_", dots$algorithm)
  
  if (model$n_channels == 1L) {
    Qs <- t(create_Q(S))
    Qm <- t(create_Q(M))
    Qd <- t(create_Q(D))
    if (need_grad) {
      objectivef <- function(pars) {
        if (any(!is.finite(exp(pars)))) {
          return(list(objective = Inf, gradient = rep(-Inf, length(pars))))
        } 
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
        list(objective = - out$loglik + 0.5 * sum(pars^2) * penalty,
             gradient = - unlist(out[-1]) + pars * penalty)
      }
    } else {
      objectivef <- function(pars) {
        if (any(!is.finite(exp(pars)))) {
          return(Inf)
        } 
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
        
        - sum(apply(out[, T_, ], 2, logSumExp)) + 0.5 * sum(pars^2) * penalty
      }
    }
  } else {
    Qs <- t(create_Q(S))
    Qm <- lapply(M, function(m) t(create_Q(m)))
    Qd <- t(create_Q(D))
    if (need_grad) {
      objectivef <- function(pars) {
        if (any(!is.finite(exp(pars)))) {
          return(list(objective = Inf, gradient = rep(-Inf, length(pars))))
        }
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
        list(objective = - out$loglik + 0.5 * sum(pars^2) * penalty,
             gradient = - unlist(out[-1]) + pars * penalty)
      }
    } else {
      objectivef <- function(pars) {
        if (any(!is.finite(exp(pars)))) {
          return(Inf)
        } 
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
        
        - sum(apply(out[, T_, ], 2, logSumExp)) + 0.5 * sum(pars^2) * penalty
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
    
    logliks <- -unlist(lapply(out, "[[", "objective"))
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
  if (model$n_channels == 1L) {
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
    loglik = -out$objective, 
    return_code = out$status,
    message = out$message,
    iterations = out$iterations,
    logliks_of_restarts = if(restarts > 0L) logliks else NULL, 
    return_codes_of_restarts = if(restarts > 0L) return_codes else NULL,
    penalty = penalty
  )
  model
}
