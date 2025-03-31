# create test data
set.seed(123)
s <- 4
n_id <- 10
n_time <- 15

data <- data.frame(
  y = factor(sample(letters[1:s], n_id * n_time, replace = TRUE)), 
  x = stats::rnorm(n_id * n_time), 
  z = stats::rnorm(n_id * n_time),
  time = rep(1:n_time, each = n_id),
  id = rep(1:n_id, n_time)
)

test_that("build_nhmm returns object of class 'nhmm'", {
  expect_error(
    model <- build_nhmm(
      "y", s, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, data = data, 
      time = "time", id = "id", state_names = 1:s,
    ),
    NA
  )
  expect_s3_class(
    model,
    "nhmm"
  )
})
test_that("estimate_nhmm returns object of class 'nhmm'", {
  expect_error(
    fit <- estimate_nhmm(
      "y", s, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, data = data, time = "time", id = "id"),
    NA
  )
  expect_s3_class(
    fit,
    "nhmm"
  )
})
test_that("estimate_nhmm errors with missing 'n_states' argument", {
  expect_error(
    estimate_nhmm("y", data = data, time = "time", id = "id"),
    "Argument `n\\_states` must be a single integer larger than 1\\."
  )
})
test_that("estimate_nhmm errors with incorrect formulas", {
  expect_error(
    estimate_nhmm("y", n_states = 3, initial_formula = 5,
                  data = data, time = "time", id = "id"),
    "Argument `initial\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_nhmm("y", n_states = 3, transition_formula = 5,
                  data = data, time = "time", id = "id"),
    "Argument `transition\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_nhmm("y", n_states = 3, emission_formula = "a",
                  data = data, time = "time", id = "id"),
    "Argument `emission\\_formula` must be a <formula> object\\."
  )
})
test_that("estimate_nhmm errors with incorrect observations", {
  expect_error(
    estimate_nhmm(list(1, "a"), s, data = data, time = "time", id = "id"),
    paste0("Argument `responses` must be a character vector defining the response variable\\(s\\) in the `data`.")
  )
})
