test_that("Gradients for singlechannel-NHMM are correct", {
  set.seed(123)
  M <- 4
  S <- 3
  n_id <- 5
  n_time <- 10
  
  data <- data.table(
    y = sample(factor(letters[1:M]), n_id * n_time, replace = TRUE), 
    x = stats::rnorm(n_id * n_time), 
    z = stats::rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  data[data$id < 3 & data$time > 6, c("y", "x", "z")] <- NA
  data[data$time > 9, c("y", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y[10:25] <- NA
  data <- stats::na.omit(data, cols = c("x", "z"))
  model <- build_nhmm(
    S, initial_formula = ~ x, transition_formula = ~ z,
    emission_formula = y ~ x, data = data, time = "time", id = "id")
  
  np_pi <- attr(model, "np_pi")
  np_A <- attr(model, "np_A")
  np_B <- attr(model, "np_B")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  obs <- create_obs(model)
  pars <- stats::rnorm(np_pi + np_A + np_B)
  
  f <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
    eta_B <- create_eta_B_nhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
    )
    -Rcpp_log_objective_nhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, TRUE, eta_pi, eta_A, eta_B
    )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
    eta_B <- create_eta_B_nhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
    )
    -unname(unlist(Rcpp_log_objective_nhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, TRUE, eta_pi, eta_A, eta_B)[-1])
    )
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

test_that("Gradients for multichannel-NHMM are correct", {
  set.seed(123)
  M <- c(2, 5)
  S <- 3
  n_id <- 5
  n_time <- 15
  data <- data.table(
    y1 = sample(factor(letters[1:M[1]]), n_id * n_time, replace = TRUE), 
    y2 = sample(factor(LETTERS[1:M[2]]), n_id * n_time, replace = TRUE), 
    x = stats::rnorm(n_id * n_time), 
    z = stats::rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  
  data[data$id < 3 & data$time > 6, c("y1", "y2", "x", "z")] <- NA
  data[data$time > 9, c("y1", "y2", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y1[10:25] <- NA
  data$y2[10:35] <- NA
  data <- stats::na.omit(data, cols = c("x", "z"))
  model <- build_nhmm(
    S, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = c(y1, y2) ~ z, data = data, time = "time", id = "id")
  
  np_pi <- attr(model, "np_pi")
  np_A <- attr(model, "np_A")
  np_B <- attr(model, "np_B")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  obs <- create_obs(model)
  pars <- stats::rnorm(np_pi + np_A + np_B)
  
  f <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
    eta_B <- create_eta_B_nhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
    )
    -Rcpp_log_objective_nhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, io(X_pi), io(X_A), io(X_B),
      iv(X_A), iv(X_B), tv(X_A), tv(X_B), eta_pi, eta_A, eta_B
    )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(np_pi)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[np_pi + seq_len(np_A)], S, K_A)
    eta_B <- create_eta_B_nhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B
    )
    -unname(unlist(Rcpp_log_objective_nhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, io(X_pi), io(X_A), io(X_B),
      iv(X_A), iv(X_B), tv(X_A), tv(X_B), eta_pi, eta_A, eta_B
    )[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

test_that("Gradients for singlechannel-MNHMM are correct", {
  set.seed(123)
  M <- 4
  S <- 3
  D <- 2
  n_id <- 5
  n_time <- 10
  data <- data.table(
    y = sample(factor(letters[1:M]), n_id * n_time, replace = TRUE), 
    x = stats::rnorm(n_id * n_time), 
    z = stats::rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  data[data$id < 3 & data$time > 6, c("y", "x", "z")] <- NA
  data[data$time > 9, c("y", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y[10:25] <- NA
  data <- stats::na.omit(data, cols = c("x", "z"))
  model <- build_mnhmm(
    S, D, emission_formula = y ~ z, initial_formula = ~ x, 
    transition_formula = ~ z, cluster_formula = ~ z, data = data, 
    time = "time", id = "id"
  )
  
  np_pi <- attr(model, "np_pi")
  np_A <- attr(model, "np_A")
  np_B <- attr(model, "np_B")
  np_omega <- attr(model, "np_omega")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  X_omega <- model$X_omega
  K_omega <- K(X_omega)
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  obs <- create_obs(model)
  pars <- stats::rnorm(np_pi + np_A + np_B + np_omega)
  f <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
    eta_B <- create_eta_B_mnhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
    )
    -Rcpp_log_objective_mnhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, X_omega, 
      io(X_pi), io(X_A), io(X_B), io(X_omega),
      iv(X_A), iv(X_B), tv(X_A), tv(X_B),
      eta_pi, eta_A, eta_B, eta_omega
    )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
    eta_B <- create_eta_B_mnhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[np_pi + np_A + np_B + seq_len(np_omega)],  D, K_omega
    )
    -unname(unlist(Rcpp_log_objective_mnhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, X_omega, 
      io(X_pi), io(X_A), io(X_B), io(X_omega),
      iv(X_A), iv(X_B), tv(X_A), tv(X_B),
      eta_pi, eta_A, eta_B, eta_omega)[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

test_that("Gradients for multichannel-MNHMM are correct", {
  set.seed(123)
  M <- c(2, 5)
  S <- 3
  D <- 4
  n_id <- 5
  n_time <- 15
  data <- data.table(
    y1 = sample(factor(letters[1:M[1]]), n_id * n_time, replace = TRUE), 
    y2 = sample(factor(LETTERS[1:M[2]]), n_id * n_time, replace = TRUE), 
    x = stats::rnorm(n_id * n_time), 
    z = stats::rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  
  data[data$id < 3 & data$time > 6, c("y1", "y2", "x", "z")] <- NA
  data[data$time > 9, c("y1", "y2", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y1[10:25] <- NA
  data$y2[10:35] <- NA
  data <- stats::na.omit(data, cols = c("x", "z"))
  model <- build_mnhmm(
    S, D, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = c(y1, y2) ~ z, cluster_formula = ~ z, data = data, 
    time = "time", id = "id")
  
  np_pi <- attr(model, "np_pi")
  np_A <- attr(model, "np_A")
  np_B <- attr(model, "np_B")
  np_omega <- attr(model, "np_omega")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  X_omega <- model$X_omega
  K_omega <- K(X_omega)
  K_pi <- K(X_pi)
  K_A <- K(X_A)
  K_B <- K(X_B)
  obs <- create_obs(model)
  pars <- stats::rnorm(np_pi + np_A + np_B + np_omega)
  f <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
    eta_B <- create_eta_B_mnhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
    )
    -Rcpp_log_objective_mnhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, X_omega, 
      io(X_pi), io(X_A), io(X_B), io(X_omega),
      iv(X_A), iv(X_B), tv(X_A), tv(X_B),
      eta_pi, eta_A, eta_B, eta_omega
    )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(np_pi)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[np_pi + seq_len(np_A)], S, K_A, D)
    eta_B <- create_eta_B_mnhmm(
      pars[np_pi + np_A + seq_len(np_B)], S, M, K_B, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[np_pi + np_A + np_B + seq_len(np_omega)], D, K_omega
    )
    -unname(unlist(Rcpp_log_objective_mnhmm(
      obs, model$sequence_lengths, M, X_pi, X_A, X_B, X_omega, 
      io(X_pi), io(X_A), io(X_B), io(X_omega),
      iv(X_A), iv(X_B), tv(X_A), tv(X_B),
      eta_pi, eta_A, eta_B, eta_omega)[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

