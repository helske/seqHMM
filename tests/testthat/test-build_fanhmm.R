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
  y = droplevels(unlist(obs)), 
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
    "For FAN-HMM, the response variable `observations` must be in the `data`."
  )
  expect_error(
    model <- build_fanhmm(
      "y", s, initial_formula = ~ x, transition_formula = ~z,
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
    "fanhmm"
  )
  expect_s3_class(
    fit,
    "nhmm"
  )
})
