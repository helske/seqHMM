# create test data
set.seed(123)
s <- 4
d <- 3
n_id <- 10
n_time <- 15
data <- data.frame(
  y = factor(sample(letters[1:s], n_id * n_time, replace = TRUE)), 
  y2 = factor(sample(letters[1:(s + 1)], n_id * n_time, replace = TRUE)),
  x = stats::rnorm(n_id * n_time), 
  z = stats::rnorm(n_id * n_time),
  time = rep(1:n_time, each = n_id),
  id = rep(1:n_id, n_time)
)
test_that("build_mnhmm returns object of class 'mnhmm'", {
  expect_error(
    model <- build_mnhmm(
      "y", s, d, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, cluster_formula = ~ x, data = data, 
      time = "time", id = "id", state_names = 1:s,
      cluster_names = 1:d
    ),
    NA
  )
  expect_s3_class(
    model,
    "mnhmm"
  )
  expect_error(
    cluster_names(model) <- seq_len(d),
    NA
  )
  expect_equal(
    cluster_names(model),
    seq_len(d)
  )
})
test_that("estimate_mnhmm returns object of class 'mnhmm'", {
  expect_error(
    fit <- estimate_mnhmm(
      "y", s, d, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, cluster_formula = ~ x,
      data = data, time = "time", id = "id", maxeval = 1,
      method = "EM"),
    NA
  )
  expect_s3_class(
    fit,
    "mnhmm"
  )
  expect_error(
    fit <- estimate_mnhmm(
      c("y", "y2"), s, d, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, cluster_formula = ~ x,
      data = data, time = "time", id = "id", maxeval = 1, method = "DNM"),
    NA
  )
  expect_s3_class(
    fit,
    "mnhmm"
  )
})
test_that("estimate_mnhmm errors with missing 'n_clusters' argument", {
  expect_error(
    estimate_mnhmm("y"),
    "Argument `n\\_clusters` must be a single positive integer larger than 1\\."
  )
})
test_that("estimate_mnhmm errors with missing 'n_states' argument", {
  expect_error(
    estimate_mnhmm("y", n_clusters = 3, data = data),
    "Argument `n\\_states` must be a single integer larger than 1\\."
  )
})
test_that("estimate_mnhmm errors with incorrect formulas", {
  expect_error(
    estimate_mnhmm("y", n_states = 3, n_clusters = 2, initial_formula = 5,
                   data = data),
    "Argument `initial\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_mnhmm("y", n_states = 3, n_clusters = 2, transition_formula = 5,
                   data = data),
    "Argument `transition\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_mnhmm("y", n_states = 3, n_clusters = 2, emission_formula = "a",
                   data = data),
    "Argument `emission\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_mnhmm("y", n_states = 3, n_clusters = 2, cluster_formula = 5,
                   data = data),
    "Argument `cluster\\_formula` must be a <formula> object\\."
  )
})
test_that("estimate_mnhmm errors with missing data arguments", {
  expect_error(
    estimate_mnhmm("y", s, d, initial_formula = ~ x, data = "a"),
    "Argument `data` must be a <data.frame> object."
  )
  expect_error(
    estimate_mnhmm("y", 3, 2, initial_formula = ~ z, transition_formula = ~x,
                   data = data),
    "Argument `time` is missing."
  )
  expect_error(
    estimate_mnhmm("y", 3, 5, emission_formula = ~x,
                   data = data, time = "time", id = 2),
    "Argument `id` must be a single character string."
  )
  expect_error(
    estimate_mnhmm("y", 3, 5, emission_formula = ~x,
                   data = data, id = "id", time = -4),
    "Argument `time` must be a single character string."
  )
})
test_that("estimate_mnhmm errors with incorrect observations", {
  expect_error(
    estimate_mnhmm(list(1,"a"), s, d, data = data, time = "time", id = "id"),
    "Argument `responses` must be a character vector defining the response variable\\(s\\) in the `data`."
  )
})

test_that("build_mnhmm works with missing observations", {
  data <- data[data$time != 5, ]
  data$y[50:55] <- NA
  expect_error(
    model <- estimate_mnhmm(
      "y", s, d, data = data, time = "time", id = "id", maxeval = 1,
      maxeval_em_dnm = 1),
    NA
  )
  expect_equal(
    which(is.na(model$data$y)),
    c(5L, 7L, 20L, 22L, 35L, 37L, 50L, 52L, 65L, 67L, 80L, 95L, 110L, 
      125L, 140L, 141L)
  )
})
