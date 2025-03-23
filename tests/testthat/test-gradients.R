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
    "y", S, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ x, data = data, time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B)
  obs <- create_obsArray(model)
  pars <- stats::rnorm(n_i + n_s + n_o)
  
  f <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
    eta_B <- create_eta_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B
    )
    -log_objective_nhmm_singlechannel(
      eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, model$sequence_lengths, 
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE
    )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
    eta_B <- create_eta_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B
    )
    -unname(unlist(log_objective_nhmm_singlechannel(
      eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, model$sequence_lengths, 
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)[-1])
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
    c("y1", "y2"), S, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ z, data = data, time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B)
  obs <- create_obsArray(model)
  pars <- stats::rnorm(n_i + n_s + n_o)
  
  f <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
    eta_B <- create_eta_multichannel_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B
    )
    -log_objective_nhmm_multichannel(
      eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, model$sequence_lengths, 
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE
    )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_pi)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_A)
    eta_B <- create_eta_multichannel_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B
    )
    -unname(unlist(log_objective_nhmm_multichannel(
      eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, model$sequence_lengths, 
      FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE
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
    "y", S, D, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ z, cluster_formula = ~ z, data = data, 
    time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  n_d <- attr(model, "np_omega")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  X_omega <- model$X_omega
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B)
  K_omega <- nrow(X_omega)
  obs <- create_obsArray(model)
  pars <- stats::rnorm(n_i + n_s + n_o + n_d)
  f <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
    eta_B <- create_eta_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
    )
    -log_objective_mnhmm_singlechannel(
      eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, 
      model$sequence_lengths, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, 
      TRUE
    )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
    eta_B <- create_eta_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)],  D, K_omega
    )
    -unname(unlist(log_objective_mnhmm_singlechannel(
      eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, 
      model$sequence_lengths, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, 
      TRUE)[-1]))
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
    c("y1", "y2"), S, D, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ z, cluster_formula = ~ z, data = data, 
    time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  n_d <- attr(model, "np_omega")
  X_pi <- model$X_pi
  X_A <- model$X_A
  X_B <- model$X_B
  X_omega <- model$X_omega
  K_pi <- nrow(X_pi)
  K_A <- nrow(X_A)
  K_B <- nrow(X_B)
  K_omega <- nrow(X_omega)
  obs <- create_obsArray(model)
  pars <- stats::rnorm(n_i + n_s + n_o + n_d)
  f <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
    eta_B <- unlist(
      create_eta_multichannel_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
      ),
      recursive = FALSE
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
    )
    -log_objective_mnhmm_multichannel(
      eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, 
      model$sequence_lengths, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, 
      TRUE
      )$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_pi, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_A, D)
    eta_B <- unlist(
      create_eta_multichannel_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_B, D
      ),
      recursive = FALSE
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_omega
    )
    -unname(unlist(log_objective_mnhmm_multichannel(
      eta_omega, X_omega, eta_pi, X_pi, eta_A, X_A, eta_B, X_B, obs, 
      model$sequence_lengths, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, 
      TRUE)[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

