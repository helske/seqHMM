# create test data
set.seed(123)
s <- 4
n_id <- 10
n_time <- 15
obs <- suppressMessages(seqdef(
  matrix(
    sample(letters[1:s], n_id * n_time, replace = TRUE), 
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

test_that("build_fanhmm returns object of class 'fanhmm'", {
  expect_error(
    model <- build_fanhmm(
      obs, s, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, autoregression_formula = ~ x, 
      feedback_formula = ~ 1, data = data, 
      time = "time", id = "id", state_names = 1:s
    ),
    NA
  )
  expect_s3_class(
    model,
    "fanhmm"
  )
})
test_that("estimate_fanhmm returns object of class 'fanhmm'", {
  expect_error(
    fit <- estimate_fanhmm(
      "y", s, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, autoregression_formula = ~ 1, 
      data = data, time = "time", id = "id", maxeval = 1),
    NA
  )
  expect_s3_class(
    fit,
    "nhmm"
  )
})
test_that("estimate_fanhmm errors with missing 'n_states' argument", {
  expect_error(
    estimate_fanhmm(obs),
    "Argument `n\\_states` must be a single integer larger than 1\\."
  )
})
test_that("estimate_fanhmm errors with incorrect formulas", {
  expect_error(
    estimate_fanhmm(obs, n_states = 3, initial_formula = 5),
    "Argument `initial\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_fanhmm(obs, n_states = 3, transition_formula = 5),
    "Argument `transition\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_fanhmm(obs, n_states = 3, emission_formula = "a"),
    "Argument `emission\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_fanhmm(obs, n_states = 3, feedback_formula = "a"),
    "Argument `feedback\\_formula` must be a <formula> object\\."
  )
})
test_that("estimate_fanhmm errors with missing data arguments", {
  expect_error(
    estimate_fanhmm(obs, s, initial_formula = ~ x),
    "Argument `data` must be a <data.frame> object."
  )
  expect_error(
    estimate_fanhmm(obs, 3, initial_formula = ~ z, transition_formula = ~x,
                  data = data),
    "Argument `time` must be a single character string."
  )
  expect_error(
    estimate_fanhmm(obs, 3, emission_formula = ~x,
                  data = data, time = "time"),
    "Argument `id` must be a single character string."
  )
  expect_error(
    estimate_fanhmm(obs, 3, emission_formula = ~x,
                  data = data, id = "id"),
    "Argument `time` must be a single character string."
  )
})

test_that("build_fanhmm works with missing observations", {
  data <- data[data$time != 5, ]
  data$y[50:55] <- NA
  expect_error(
    model <- estimate_fanhmm(
      "y", s, data = data, time = "time", id = "id", lambda = 1),
    "FAN-HMM does not support missing values in the observations."
  )
})
