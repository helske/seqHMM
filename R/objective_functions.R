make_objective_nhmm <- function(model, lambda = 0, need_grad = TRUE) {
  M <- model$n_symbols
  S <- model$n_states
  T_ <- model$length_of_sequences
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
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  obs <- create_obs(model)
  
  W_X_B <- model$W_X_B
  prior_obs <- model$prior_obs
  use_fanhmm <- inherits(model, "fanhmm") && !identical(prior_obs, 0L)
  
  if (need_grad) {
    function(pars) {
      
      eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
      eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
      eta_B <- create_eta_B_nhmm(pars[np_pi + np_A + seq_len(np_B)], S, M, K_B)
      
      if (use_fanhmm) {
        out <- Rcpp_log_objective_fanhmm(
          obs, Ti, M, X_pi, X_A, X_B,
          icpt_only_pi, icpt_only_A, icpt_only_B,
          iv_A, iv_B, tv_A, tv_B,
          eta_pi, eta_A, eta_B, prior_obs, W_X_B
        )
      } else {
        out <- Rcpp_log_objective_nhmm(
          obs, Ti, M, X_pi, X_A, X_B,
          icpt_only_pi, icpt_only_A, icpt_only_B,
          iv_A, iv_B, tv_A, tv_B,
          eta_pi, eta_A, eta_B
        )
      }
      
      list(
        objective = -(out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs,
        gradient = -(unlist(out[-1]) - lambda * pars) / n_obs
      )
    }
  } else {
    function(pars) {
      
      eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
      eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
      eta_B <- create_eta_B_nhmm(pars[np_pi + np_A + seq_len(np_B)], S, M, K_B)
      
      if (use_fanhmm) {
        ll <- Rcpp_loglik_fanhmm(
          obs, Ti, M, X_pi, X_A, X_B,
          icpt_only_pi, icpt_only_A, icpt_only_B,
          iv_A, iv_B, tv_A, tv_B,
          eta_pi, eta_A, eta_B, prior_obs, W_X_B
        )
      } else {
        ll <- Rcpp_loglik_nhmm(
          obs, Ti, M, X_pi, X_A, X_B,
          icpt_only_pi, icpt_only_A, icpt_only_B,
          iv_A, iv_B, tv_A, tv_B,
          eta_pi, eta_A, eta_B
        )
      }
      -(sum(ll) - 0.5 * lambda * sum(pars^2)) / n_obs
    }
  }
}
make_objective_mnhmm <- function(model, lambda = 0, need_grad = TRUE) {
  M <- model$n_symbols
  S <- model$n_states
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
  K_omega <- K(X_omega)
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  Ti <- model$sequence_lengths
  n_obs <- nobs(model)
  obs <- create_obs(model)
  
  W_X_B <- model$W_X_B
  prior_obs <- model$prior_obs
  use_fanhmm <- inherits(model, "fanhmm") && !identical(prior_obs, 0L)
  
  if (need_grad) {
    function(pars) {
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
      if (use_fanhmm) {
        out <- Rcpp_log_objective_mfanhmm(
          obs, Ti, M, X_pi, X_A, X_B, X_omega,
          icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
          iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega,
          prior_obs, W_X_B
        )
      } else {
        out <- Rcpp_log_objective_mnhmm(
          obs, Ti, M, X_pi, X_A, X_B, X_omega,
          icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
          iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega
        )
      }
      list(
        objective = -(out$loglik - 0.5 * lambda * sum(pars^2)) / n_obs,
        gradient = -(unlist(out[-1]) - lambda * pars) / n_obs
      )
    }
  } else {
    function(pars) {
      
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
      if (use_fanhmm) {
        ll <- Rcpp_loglik_mfanhmm(
          obs, Ti, M, X_pi, X_A, X_B, X_omega,
          icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
          iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega,
          prior_obs, W_X_B
        )
      } else {
        ll <- Rcpp_loglik_mnhmm(
          obs, Ti, M, X_pi, X_A, X_B, X_omega,
          icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega, 
          iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega
        )
      }
      
      - (sum(ll) - 0.5 * lambda * sum(pars^2)) / n_obs
    }
  }
}
