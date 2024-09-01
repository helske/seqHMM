# create test data
set.seed(123)
s <- 4
d <- 3
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
test_that("build_mnhmm returns object of class 'mnhmm'", {
  expect_error(
    model <- build_mnhmm(
      obs, s, d, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, cluster_formula = ~ x, data = data, 
      time = "time", id = "id", state_names = 1:s, channel_names = "obs", 
      cluster_names = 1:d
    ),
    NA
  )
  expect_s3_class(
    model,
    "mnhmm"
  )
})
test_that("estimate_mnhmm returns object of class 'mnhmm'", {
  expect_error(
    fit <- estimate_mnhmm(
      "y", s, d, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, cluster_formula = ~ x,
      data = data, time = "time", id = "id", 
      iter = 0),
    NA
  )
  expect_s3_class(
    fit,
    "mnhmm"
  )
})
test_that("estimate_mnhmm errors with missing 'n_clusters' argument", {
  expect_error(
    estimate_mnhmm(obs),
    "Argument `n\\_clusters` must be a single positive integer larger than 1\\."
  )
})
test_that("estimate_mnhmm errors with missing 'n_states' argument", {
  expect_error(
    estimate_mnhmm(obs, n_clusters = 3),
    "Argument `n\\_states` must be a single positive integer\\."
  )
})
test_that("estimate_mnhmm errors with incorrect formulas", {
  expect_error(
    estimate_mnhmm(obs, n_states = 3, n_clusters = 2, initial_formula = 5),
    "Argument `initial\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_mnhmm(obs, n_states = 3, n_clusters = 2, transition_formula = 5),
    "Argument `transition\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_mnhmm(obs, n_states = 3, n_clusters = 2, emission_formula = "a"),
    "Argument `emission\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_mnhmm(obs, n_states = 3, n_clusters = 2, cluster_formula = 5),
    "Argument `cluster\\_formula` must be a <formula> object\\."
  )
})
test_that("estimate_mnhmm errors with missing data arguments", {
  expect_error(
    estimate_mnhmm(obs, s, d, initial_formula = ~ x),
    "Argument `data` must be a <data.frame> object."
  )
  expect_error(
    estimate_mnhmm(obs, 3, 2, initial_formula = ~ z, transition_formula = ~x,
                   data = data),
    "Argument `time` must be a single character string."
  )
  expect_error(
    estimate_mnhmm(obs, 3, 5, emission_formula = ~x,
                   data = data, time = "time"),
    "Argument `id` must be a single character string."
  )
  expect_error(
    estimate_mnhmm(obs, 3, 5, emission_formula = ~x,
                   data = data, id = "id"),
    "Argument `time` must be a single character string."
  )
})
test_that("estimate_mnhmm errors with incorrect observations", {
  expect_error(
    estimate_mnhmm(list(1,"a"), s, d),
    paste0("`observations` should be a <stslist> object created with ",
           "`seqdef\\(\\)`, a <list> of <stslist> objects, or a <character> ",
           "vector containing names of the response variables in `data`."
    )
  )
})
test_that("build_mnhmm works with vector of characters as observations", {
  expect_error(
    estimate_mnhmm("y", s, d, data = data, time = "time", id = "id", iter = 0),
    NA
  )
})
