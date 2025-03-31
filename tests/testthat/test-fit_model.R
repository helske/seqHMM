
test_that("'fit_model' works for 'hmm'", {
  set.seed(123)
  emission_probs <- matrix(c(0.5, 0.2, 0.5, 0.8), 2, 2)
  transition_probs <- matrix(c(5 / 6, 1 / 6, 1 / 6, 5 / 6), 2, 2)
  initial_probs <- c(1, 0)
  expect_error(
    sim <- simulate_hmm(
      n_sequences = 10, initial_probs = initial_probs,
      transition_probs = transition_probs,
      emission_probs = emission_probs,
      sequence_length = 20
    ),
    NA
  )
  
  model <- build_hmm(
    observations = sim$observations,
    transition_probs = transition_probs,
    emission_probs = emission_probs,
    initial_probs = initial_probs
  )
  set.seed(1)
  expect_error(
    fit2 <- fit_model(
      model, em_step = TRUE, global_step = TRUE,
      local_step = TRUE,
      control_em = list(restart = list(times = 2), maxeval = 10),
      control_global = list(maxeval = 10),
      control_local = list(maxeval = 10)
    ),
    NA
  )
  set.seed(1)
  expect_error(
    fit1 <- fit_model(
      model, em_step = TRUE, global_step = TRUE,
      local_step = TRUE,
      control_em = list(restart = list(times = 2), maxeval = 10),
      control_global = list(maxeval = 10),
      control_local = list(maxeval = 10),
      log_space = FALSE
    ),
    NA
  )
  expect_equal(fit1$model, fit2$model)
  expect_equal(logLik(fit1$model), logLik(fit2$model, log_space = FALSE))
  expect_equal(
    logLik(fit1$model, threads = 2L)[1], 
    sum(logLik(fit1$model, partials = TRUE))
  )
})

test_that("'fit_model' works for 'mhmm'", {
  set.seed(123)
  emission_probs_1 <- matrix(c(0.75, 0.05, 0.25, 0.95), 2, 2)
  emission_probs_2 <- matrix(c(0.1, 0.8, 0.9, 0.2), 2, 2)
  colnames(emission_probs_1) <- colnames(emission_probs_2) <-
    c("heads", "tails")
  
  transition_probs_1 <- matrix(c(9, 0.1, 1, 9.9) / 10, 2, 2)
  transition_probs_2 <- matrix(c(35, 1, 1, 35) / 36, 2, 2)
  rownames(emission_probs_1) <- rownames(transition_probs_1) <-
    colnames(transition_probs_1) <- c("coin 1", "coin 2")
  rownames(emission_probs_2) <- rownames(transition_probs_2) <-
    colnames(transition_probs_2) <- c("coin 3", "coin 4")
  
  initial_probs_1 <- c(1, 0)
  initial_probs_2 <- c(1, 0)
  
  n <- 30
  covariate_1 <- runif(n)
  covariate_2 <- sample(c("A", "B"),
                        size = n, replace = TRUE,
                        prob = c(0.3, 0.7)
  )
  dataf <- data.frame(covariate_1, covariate_2)
  
  coefs <- cbind(cluster_1 = c(0, 0, 0), cluster_2 = c(-1.5, 3, -0.7))
  rownames(coefs) <- c("(Intercept)", "covariate_1", "covariate_2B")
  
  expect_error(
    sim <- simulate_mhmm(
      n = n, initial_probs = list(initial_probs_1, initial_probs_2),
      transition_probs = list(transition_probs_1, transition_probs_2),
      emission_probs = list(emission_probs_1, emission_probs_2),
      sequence_length = 20, formula = ~ covariate_1 + covariate_2,
      data = dataf, coefficients = coefs
    ),
    NA
  )
  
  expect_error(
    model <- build_mhmm(
      sim$observations,
      initial_probs = list(initial_probs_1, initial_probs_2),
      transition_probs = list(transition_probs_1, transition_probs_2),
      emission_probs = list(emission_probs_1, emission_probs_2),
      formula = ~ covariate_1 + covariate_2,
      data = dataf),
    NA
  )
  
  set.seed(1)
  expect_error(
    fit1 <- fit_model(
      model, em_step = TRUE, global_step = TRUE,
      local_step = TRUE,
      control_em = list(restart = list(times = 2), maxeval = 10),
      control_global = list(maxeval = 10),
      control_local = list(maxeval = 10)
    ),
    NA
  )
  
  set.seed(1)
  expect_error(
    fit2 <- fit_model(
      model, em_step = TRUE, global_step = TRUE,
      local_step = TRUE,
      control_em = list(restart = list(times = 2), maxeval = 10),
      control_global = list(maxeval = 10),
      control_local = list(maxeval = 10),
      log_space = FALSE
    ),
    NA
  )
  expect_equal(fit1$model, fit2$model)
  expect_equal(logLik(fit1$model), logLik(fit2$model, log_space = FALSE))
  expect_equal(
    logLik(fit1$model, threads = 2L)[1], 
    sum(logLik(fit1$model, partials = TRUE))
  )
})
