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

test_that("estimate_nhmm returns object of class 'fanhmm'", {
  expect_error(
    model <- estimate_nhmm(
      s, emission_formula =  ~ z * lag(y),
      initial_formula = ~ x, transition_formula = ~z + lag(y),
      data = data, time = "time", id = "id", state_names = 1:s
    ),
    paste0(
      "`emission_formula` must contain the response variable\\(s\\) on the ",
      "left-hand side of the <formula> object\\(s\\)."
    )
  )
  expect_error(
    model <- build_nhmm(
      s, emission_formula = y ~ z * lag(y),
      initial_formula = ~ x, transition_formula = ~z + lag(y),
      data = data, time = "time", id = "id", state_names = 1:s
    ),
    NA
  )
  expect_s3_class(
    model,
    c("fanhmm", "nhmm")
  )
  expect_error(
    fit <- estimate_nhmm(
      s, emission_formula = y ~ z * lag(y),
      initial_formula = ~ x, transition_formula = ~z + lag(y),
      data = data, time = "time", id = "id", maxeval = 1
    ),
    NA
  )
  expect_s3_class(
    model,
    c("fanhmm", "nhmm")
  )
  expect_error(
    model <- build_mnhmm(
      s, 2, emission_formula = y ~ z * lag(y),
      initial_formula = ~ x, transition_formula = ~z + lag(y),
      cluster_formula = ~ z,
      data = data, time = "time", id = "id", state_names = 1:s
    ),
    NA
  )
  expect_s3_class(
    model,
    c("fanhmm", "mnhmm")
  )
})

