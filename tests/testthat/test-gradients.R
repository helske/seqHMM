test_that("Gradients for singlechannel-NHMM are correct", {
  set.seed(123)
  M <- 4
  S <- 3
  n_id <- 5
  n_time <- 10
  obs <- suppressMessages(seqdef(
    matrix(
      sample(letters[1:M], n_id * n_time, replace = TRUE), 
      n_id, n_time
    )
  ))
  data <- data.frame(
    y = unlist(obs), 
    x = rnorm(n_id * n_time), 
    z = rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  data[data$id < 3 & data$time > 6, c("y", "x", "z")] <- NA
  data[data$time > 9, c("y", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y[10:25] <- NA
  model <- build_nhmm(
    "y", S, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ x, data = data, time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  X_i <- model$X_initial
  X_s <- model$X_transition
  X_o <- model$X_emission
  K_i <- nrow(X_i)
  K_s <- nrow(X_s)
  K_o <- nrow(X_o)
  obs <- create_obsArray(model)
  obs <- array(obs, dim(obs)[2:3])
  pars <- rnorm(n_i + n_s + n_o)
 
  Qs <- t(create_Q(S))
  Qm <- t(create_Q(M))
  f <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_i)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    eta_B <- create_eta_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -log_objective_nhmm_singlechannel(
      Qs, Qm,
      eta_pi, X_i,
      eta_A, X_s,
      eta_B, X_o,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_i)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    eta_B <- create_eta_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -unname(unlist(log_objective_nhmm_singlechannel(
      Qs, Qm, eta_pi, X_i, eta_A, X_s, eta_B, X_o,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

test_that("Gradients for multichannel-NHMM are correct", {
  set.seed(123)
  M <- c(2, 5)
  S <- 3
  n_id <- 5
  n_time <- 15
  obs1 <- suppressMessages(seqdef(
    matrix(
      sample(letters[1:M[1]], n_id * n_time, replace = TRUE), 
      n_id, n_time
    )
  ))
  obs2 <- suppressMessages(seqdef(
    matrix(
      sample(LETTERS[1:M[2]], n_id * n_time, replace = TRUE), 
      n_id, n_time
    )
  ))
  data <- data.frame(
    y1 = unlist(obs1), 
    y2 = unlist(obs2), 
    x = rnorm(n_id * n_time), 
    z = rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  
  data[data$id < 3 & data$time > 6, c("y1", "y2", "x", "z")] <- NA
  data[data$time > 9, c("y1", "y2", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y1[10:25] <- NA
  data$y2[10:35] <- NA
  
  model <- build_nhmm(
    c("y1", "y2"), S, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ z, data = data, time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  X_i <- model$X_initial
  X_s <- model$X_transition
  X_o <- model$X_emission
  K_i <- nrow(X_i)
  K_s <- nrow(X_s)
  K_o <- nrow(X_o)
  obs <- create_obsArray(model)
  pars <- rnorm(n_i + n_s + n_o)
  Qs <- t(create_Q(S))
  Qm <- lapply(M, function(m) t(create_Q(m)))
  
  f <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_i)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    eta_B <- create_eta_multichannel_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -log_objective_nhmm_multichannel(
      Qs, Qm, eta_pi, X_i, eta_A, X_s, eta_B, X_o,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_nhmm(pars[seq_len(n_i)], S, K_i)
    eta_A <- create_eta_A_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    eta_B <- create_eta_multichannel_B_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -unname(unlist(log_objective_nhmm_multichannel(
      Qs, Qm, eta_pi, X_i, eta_A, X_s, eta_B, X_o,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)[-1]))
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
  obs <- suppressMessages(seqdef(
    matrix(
      sample(letters[1:M], n_id * n_time, replace = TRUE), 
      n_id, n_time
    )
  ))
  data <- data.frame(
    y = unlist(obs), 
    x = rnorm(n_id * n_time), 
    z = rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  data[data$id < 3 & data$time > 6, c("y", "x", "z")] <- NA
  data[data$time > 9, c("y", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y[10:25] <- NA
  model <- build_mnhmm(
    "y", S, D, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ z, cluster_formula = ~ z, data = data, 
    time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  n_d <- attr(model, "np_omega")
  X_i <- model$X_initial
  X_s <- model$X_transition
  X_o <- model$X_emission
  X_d <- model$X_cluster
  K_i <- nrow(X_i)
  K_s <- nrow(X_s)
  K_o <- nrow(X_o)
  K_d <- nrow(X_d)
  obs <- create_obsArray(model)
  obs <- array(obs, dim(obs)[2:3])
  pars <- rnorm(n_i + n_s + n_o + n_d)
  Qs <- t(create_Q(S))
  Qm <- t(create_Q(M))
  Qd <- t(create_Q(D))
  f <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    eta_B <- create_eta_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    -log_objective_mnhmm_singlechannel(
      Qs, Qm, Qd, eta_pi, X_i, eta_A, X_s, eta_B, X_o, eta_omega, X_d,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    eta_B <- create_eta_B_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)],  D, K_d
    )
    -unname(unlist(log_objective_mnhmm_singlechannel(
      Qs, Qm, Qd, eta_pi, X_i, eta_A, X_s, eta_B, X_o, eta_omega, X_d,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)[-1]))
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
  obs1 <- suppressMessages(seqdef(
    matrix(
      sample(letters[1:M[1]], n_id * n_time, replace = TRUE), 
      n_id, n_time
    )
  ))
  obs2 <- suppressMessages(seqdef(
    matrix(
      sample(LETTERS[1:M[2]], n_id * n_time, replace = TRUE), 
      n_id, n_time
    )
  ))
  data <- data.frame(
    y1 = unlist(obs1), 
    y2 = unlist(obs2), 
    x = rnorm(n_id * n_time), 
    z = rnorm(n_id * n_time),
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  
  data[data$id < 3 & data$time > 6, c("y1", "y2", "x", "z")] <- NA
  data[data$time > 9, c("y1", "y2", "x", "z")] <- NA
  data$x[12:15] <- 0
  data$y1[10:25] <- NA
  data$y2[10:35] <- NA
  
  model <- build_mnhmm(
    c("y1", "y2"), S, D, initial_formula = ~ x, transition_formula = ~z,
    emission_formula = ~ z, cluster_formula = ~ z, data = data, 
    time = "time", id = "id")
  
  n_i <- attr(model, "np_pi")
  n_s <- attr(model, "np_A")
  n_o <- attr(model, "np_B")
  n_d <- attr(model, "np_omega")
  X_i <- model$X_initial
  X_s <- model$X_transition
  X_o <- model$X_emission
  X_d <- model$X_cluster
  K_i <- nrow(X_i)
  K_s <- nrow(X_s)
  K_o <- nrow(X_o)
  K_d <- nrow(X_d)
  obs <- create_obsArray(model)
  pars <- rnorm(n_i + n_s + n_o + n_d)
  Qs <- t(create_Q(S))
  Qm <- lapply(M, function(m) t(create_Q(m)))
  Qd <- t(create_Q(D))
  f <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    eta_B <- unlist(
      create_eta_multichannel_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      ),
      recursive = FALSE
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    -log_objective_mnhmm_multichannel(
      Qs, Qm, Qd, eta_pi, X_i, eta_A, X_s, eta_B, X_o, eta_omega, X_d,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)$loglik
    
  }
  g <- function(pars) {
    eta_pi <- create_eta_pi_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    eta_A <- create_eta_A_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    eta_B <- unlist(
      create_eta_multichannel_B_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      ),
      recursive = FALSE
    )
    eta_omega <- create_eta_omega_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    -unname(unlist(log_objective_mnhmm_multichannel(
      Qs, Qm, Qd, eta_pi, X_i, eta_A, X_s, eta_B, X_o, eta_omega, X_d,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, model$sequence_lengths)[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

