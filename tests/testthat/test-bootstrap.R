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
    fit <- estimate_nhmm(3, y ~ 1, data = d, time = "time", id = "id"),
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
      method = "DNM"
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
