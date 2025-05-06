test_that("boostrap works for `nhmm`", {
  set.seed(123)
  s <- 2
  n_id <- 10
  n_time <- 15
  
  d <- data.frame(
    y = factor(sample(letters[1:s], n_id * n_time, replace = TRUE)), 
    x = stats::rnorm(n_id * n_time), 
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  expect_error(
    fit <- estimate_nhmm(3, y ~ 1, data = d, time = "time", id = "id", 
                         maxeval = 2),
    NA
  )
  expect_error(
    fit <- bootstrap_coefs(fit, type = "nonparametric", nsim = 2),
    NA
  )
  expect_error(
    fit <- bootstrap_coefs(fit, type = "parametric", nsim = 2, append = TRUE),
    NA
  )
  expect_equal(length(fit$boot), 4L)
  expect_error(
    bootstrap_coefs(fit), 
    "Argument `nsim` must be a single positive integer."
  )
  expect_error(
    bootstrap_coefs(fit, "a"), 
    "Argument `nsim` must be a single positive integer."
  )
  expect_error(
    bootstrap_coefs(fit, 1, "a"), 
    'Argument `type` must be either \"nonparametric\" or \"parametric\"\\.'
  )
})

test_that("boostrap works for `mnhmm`", {
  set.seed(123)
  s <- 4
  n_id <- 30
  n_time <- 10
  
  d <- data.frame(
    y = factor(sample(letters[1:s], n_id * n_time, replace = TRUE)), 
    x = stats::rnorm(n_id * n_time), 
    time = rep(1:n_time, each = n_id),
    id = rep(1:n_id, n_time)
  )
  expect_error(
    fit <- estimate_mnhmm(
      3, 2, emission_formula = y ~ 1, data = d, time = "time", id = "id", 
      method = "DNM", maxeval = 10
    ),
    NA
  )
  expect_error(
    fit <- bootstrap_coefs(fit, type = "nonparametric", nsim = 2),
    NA
  )
  expect_error(
    fit <- bootstrap_coefs(fit, type = "parametric", nsim = 2, append = TRUE),
    NA
  )
  expect_equal(length(fit$boot), 5L)
  expect_error(
    bootstrap_coefs(fit), 
    "Argument `nsim` must be a single positive integer."
  )
  expect_error(
    bootstrap_coefs(fit, "a"), 
    "Argument `nsim` must be a single positive integer."
  )
  expect_error(
    bootstrap_coefs(fit, 1, "a"), 
    'Argument `type` must be either \"nonparametric\" or \"parametric\"\\.'
  )
})


test_that("boostrap works for `mfanhmm`", {
  set.seed(123)
  n_seq <- 100
  t_seq <- 20
  pi <- list(c(0.6, 0.4), c(0.8, 0.2))
  omega = c(0.3, 0.7)
  A <- list(diag(0.8, 2) + 0.1, matrix(c(0.95, 0.2, 0.05, 0.8), 2, 2))
  B <- list(
    list(rbind(c(0.5,0.2,0.2,0.1), c(0.1, 0.1, 0.7, 0.1))), 
    list(rbind(c(0.1,0.1,0.4,0.4), c(0.2, 0.3, 0.1, 0.4)))
  )
  expect_error(
    sim <- simulate_mnhmm(
      2, 2, y ~ lag(y),
      data = data.frame(
        id = rep(1:n_seq, each = t_seq),
        time = rep_len(1:t_seq, n_seq),
        y = factor(sample(letters[1:4], t_seq * n_seq, replace = TRUE), 
                   levels = letters[1:4])),
      id = "id", time = "time",
      coefs = list(cluster_probs = omega, initial_probs = pi, 
                   transition_probs = A, emission_probs = B)
    ),
    NA
  )
  y_prior <- prop.table(table(sim$data$y[sim$data$time == 1]))
  expect_error(
    fit <- estimate_mnhmm(
      n_states = 2, n_clusters = 2,
      data = sim$data, time = "time", id = "id", 
      emission_formula = y ~ lag(y), 
      inits = sim$model$etas, init_sd = 0, lambda = 0.1,
      method = "DNM", prior_obs = list(y = y_prior)
    ),
    NA
  )
  expect_error(
    fit <- bootstrap_coefs(fit, type = "nonparametric", nsim = 2),
    NA
  )
  expect_equal(length(fit$boot), 5L)
  expect_equal(length(fit$boot[[1]]), 2L)
})
