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
  data <- data[-10L, ]
  data$y[10:15] <- NA
  data$x[12] <- NA
  model <- build_nhmm(
    "y", S, initial_formula = ~ x, transition_formula = ~z,
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
  obs <- array(obs, dim(obs)[2:3])
  pars <- rnorm(n_i + n_s + n_o)
  
  f <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_nhmm(pars[seq_len(n_i)], S, K_i)
    gamma_A_raw <- create_gamma_A_raw_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    gamma_B_raw <- create_gamma_B_raw_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -log_objective_nhmm_singlechannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE)$loglik
    
  }
  g <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_nhmm(pars[seq_len(n_i)], S, K_i)
    gamma_A_raw <- create_gamma_A_raw_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    gamma_B_raw <- create_gamma_B_raw_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -unname(unlist(log_objective_nhmm_singlechannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE)[-1]))
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
  
  data <- data[-10L, ]
  data$y1[10:15] <- NA
  data$y2[12:20] <- NA
  data$x[12] <- NA
  
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
  
  f <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_nhmm(pars[seq_len(n_i)], S, K_i)
    gamma_A_raw <- create_gamma_A_raw_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    gamma_B_raw <- create_gamma_multichannel_B_raw_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -log_objective_nhmm_multichannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE)$loglik
    
  }
  g <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_nhmm(pars[seq_len(n_i)], S, K_i)
    gamma_A_raw <- create_gamma_A_raw_nhmm(pars[n_i + seq_len(n_s)], S, K_s)
    gamma_B_raw <- create_gamma_multichannel_B_raw_nhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o
    )
    -unname(unlist(log_objective_nhmm_multichannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE)[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

test_that("Gradients for singlechannel-NHMM are correct", {
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
  data <- data[-10L, ]
  data$y[10:15] <- NA
  data$x[12] <- NA
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
  
  f <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    gamma_A_raw <- create_gamma_A_raw_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    gamma_B_raw <- create_gamma_B_raw_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
    gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    -log_objective_mnhmm_singlechannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      gamma_omega_raw, X_d,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)$loglik
    
  }
  g <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    gamma_A_raw <- create_gamma_A_raw_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    gamma_B_raw <- create_gamma_B_raw_mnhmm(
      pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
    )
    gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)],  D, K_d
    )
    -unname(unlist(log_objective_mnhmm_singlechannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      gamma_omega_raw, X_d,
      obs, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)[-1]))
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
  
  data <- data[-10L, ]
  data$y1[10:15] <- NA
  data$x[12] <- NA
  
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
  
  f <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    gamma_A_raw <- create_gamma_A_raw_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    gamma_B_raw <- unlist(
      create_gamma_multichannel_B_raw_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      ),
      recursive = FALSE
    )
    gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    -log_objective_mnhmm_multichannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      gamma_omega_raw, X_d,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)$loglik
    
  }
  g <- function(pars) {
    gamma_pi_raw <- create_gamma_pi_raw_mnhmm(pars[seq_len(n_i)], S, K_i, D)
    gamma_A_raw <- create_gamma_A_raw_mnhmm(pars[n_i + seq_len(n_s)], S, K_s, D)
    gamma_B_raw <- unlist(
      create_gamma_multichannel_B_raw_mnhmm(
        pars[n_i + n_s + seq_len(n_o)], S, M, K_o, D
      ),
      recursive = FALSE
    )
    gamma_omega_raw <- create_gamma_omega_raw_mnhmm(
      pars[n_i + n_s + n_o + seq_len(n_d)], D, K_d
    )
    -unname(unlist(log_objective_mnhmm_multichannel(
      gamma_pi_raw, X_i,
      gamma_A_raw, X_s,
      gamma_B_raw, X_o,
      gamma_omega_raw, X_d,
      obs, M, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)[-1]))
  }
  expect_equal(g(pars), numDeriv::grad(f, pars))
})

