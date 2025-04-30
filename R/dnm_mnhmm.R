dnm_mnhmm <- function(model, inits, init_sd, restarts, lambda, bound, control, 
                      control_restart, save_all_solutions) {
  
  all_solutions <- NULL
  need_grad <- grepl("NLOPT_LD_", control$algorithm)
  objectivef <- make_objective_mnhmm(
    model, lambda, need_grad
  )
  if (restarts > 0L) {
    .fun <- function(base_init, u) {
      init <- base_init + u
      fit <- nloptr(
        x0 = init, eval_f = objectivef, lb = -rep(bound, length(init)), 
        ub = rep(bound, length(init)), opts = control_restart
      )
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
    fit <- future.apply::future_lapply(
      seq_len(restarts), \(i) .fun(base_init, u[, i]), future.seed = TRUE
    )
    logliks <- -unlist(lapply(fit, "[[", "objective")) * nobs(model)
    return_codes <- unlist(lapply(fit, "[[", "status"))
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
      init <- fit[[optimum]]$solution
    }
    if (save_all_solutions) {
      all_solutions <- fit
    }
  } else {
    init <- unlist(create_initial_values(inits, model, init_sd))
  }
  
  fit <- nloptr(
    x0 = init, eval_f = objectivef, lb = -rep(bound, length(init)), 
    ub = rep(bound, length(init)), opts = control
  )
  if (fit$status < 0) {
    warning_(
      paste("Optimization terminated due to error:", return_msg(fit$status))
    )
    loglik <- NaN
  } else {
    loglik <- -fit$objective * nobs(model)
  }
  pars <- fit$solution
  
  M <- model$n_symbols
  S <- model$n_states
  D <- model$n_clusters
  np_pi <- attr(model, "np_pi")
  np_A <- attr(model, "np_A")
  np_B <- attr(model, "np_B")
  np_omega <- attr(model, "np_omega")
  K_omega <- K(model$X_omega)
  K_pi <- K(model$X_pi)
  K_A <- K(model$X_A)
  K_B <- K(model$X_B)
  
  model$etas$eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
  model$gammas$gamma_pi <- drop(eta_to_gamma_mat_field(model$etas$eta_pi))
  model$etas$eta_A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
  model$gammas$gamma_A <- drop(eta_to_gamma_cube_field(model$etas$eta_A))
  model$etas$eta_B <- create_eta_B_mnhmm(
    pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
  )
  model$gammas$gamma_B <- split(
    eta_to_gamma_cube_2d_field(model$etas$eta_B), seq_len(D)
  )
  model$etas$eta_omega <- create_eta_omega_mnhmm(
    pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
  )
  model$gammas$gamma_omega <- eta_to_gamma_mat(model$etas$eta_omega)
  model$estimation_results <- list(
    loglik = loglik, 
    return_code = fit$status,
    message = fit$message,
    iterations = fit$iterations,
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
