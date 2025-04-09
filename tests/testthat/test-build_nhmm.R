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
      s, y ~ z, ~ x, ~z, data = data, 
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
      s, emission_formula = y ~ z, initial_formula = ~ x, 
      transition_formula = ~ z, data = data, time = "time", id = "id"),
    NA
  )
  expect_s3_class(
    fit,
    "nhmm"
  )
})
test_that("estimate_nhmm errors with missing 'n_states' argument", {
  expect_error(
    estimate_nhmm(emission_formula = y ~ 1, data = data, time = "t", id = "i"),
    "Argument `n\\_states` must be a single integer larger than 1\\."
  )
})
test_that("estimate_nhmm errors with incorrect formulas", {
  expect_error(
    estimate_nhmm(n_states = 3, data = data, time = "time", id = "id"),
    "Argument `emission_formula` is missing."
  )
  expect_error(
    estimate_nhmm(n_states = 3, emission_formula = 5, initial_formula = 5,
                  data = data, time = "time", id = "id"),
    paste0(
      "`emission_formula` must be a <formula> object or a list ",
      "of <formula> objects."
    )
  )
  expect_error(
    estimate_nhmm(n_states = 3, emission_formula =  ~ x, initial_formula = 5,
                  data = data, time = "time", id = "id"),
   paste0(
     "`emission\\_formula` must contain the response variable\\(s\\) on the ",
     "left-hand side of the <formula> object\\(s\\)."
   )
  )
  expect_error(
    estimate_nhmm(n_states = 3, emission_formula = y ~ x, initial_formula = 5,
                  data = data, time = "time", id = "id"),
    "Argument `initial_formula` must be a <formula> object."
  )
  expect_error(
    estimate_nhmm(2, emission_formula = y ~ x, transition_formula = 5,
                  data = data, time = "time", id = "id"),
    "Argument `transition_formula` must be a <formula> object."
  )
})
test_that("estimate_nhmm errors with incorrect observations", {
  expect_error(
    fit <- estimate_nhmm(
      s, emission_formula = o ~ z, initial_formula = ~ x, 
      transition_formula = ~ z, data = data, time = "time", id = "id"),
    "Can't find response variable `o` in `data`."
  )
})
