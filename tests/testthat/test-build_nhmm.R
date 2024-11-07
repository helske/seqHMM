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

test_that("build_nhmm returns object of class 'nhmm'", {
  expect_error(
    model <- build_nhmm(
      obs, s, initial_formula = ~ x, transition_formula = ~z,
      emission_formula = ~ z, data = data, 
      time = "time", id = "id", state_names = 1:s, channel_names = "obs"
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
      emission_formula = ~ z, data = data, time = "time", id = "id",
      lambda = 1),
    NA
  )
  expect_s3_class(
    fit,
    "nhmm"
  )
})
test_that("estimate_nhmm errors with missing 'n_states' argument", {
  expect_error(
    estimate_nhmm(obs),
    "Argument `n\\_states` must be a single positive integer\\."
  )
})
test_that("estimate_nhmm errors with incorrect formulas", {
  expect_error(
    estimate_nhmm(obs, n_states = 3, initial_formula = 5),
    "Argument `initial\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_nhmm(obs, n_states = 3, transition_formula = 5),
    "Argument `transition\\_formula` must be a <formula> object\\."
  )
  expect_error(
    estimate_nhmm(obs, n_states = 3, emission_formula = "a"),
    "Argument `emission\\_formula` must be a <formula> object\\."
  )
})
test_that("estimate_nhmm errors with missing data arguments", {
  expect_error(
    estimate_nhmm(obs, s, initial_formula = ~ x),
    "Argument `data` must be a <data.frame> object."
  )
  expect_error(
    estimate_nhmm(obs, 3, initial_formula = ~ z, transition_formula = ~x,
                  data = data),
    "Argument `time` must be a single character string."
  )
  expect_error(
    estimate_nhmm(obs, 3, emission_formula = ~x,
                  data = data, time = "time"),
    "Argument `id` must be a single character string."
  )
  expect_error(
    estimate_nhmm(obs, 3, emission_formula = ~x,
                  data = data, id = "id"),
    "Argument `time` must be a single character string."
  )
})
test_that("estimate_nhmm errors with incorrect observations", {
  expect_error(
    estimate_nhmm(list(1, "a"), s),
    paste0("`observations` should be a <stslist> object created with ",
           "`seqdef\\(\\)`, a <list> of <stslist> objects, or a <character> ",
           "vector containing names of the response variables in `data`."
    )
  )
})
test_that("build_nhmm works with vector of characters as observations", {
  expect_error(
    estimate_nhmm("y", s, data = data, time = "time", id = "id", lambda = 1),
    NA
  )
})
test_that("build_nhmm works with missing observations", {
  data <- data[data$time != 5, ]
  data$y[50:55] <- NA
  expect_error(
    model <- estimate_nhmm(
      "y", s, data = data, time = "time", id = "id", lambda = 1),
    NA
  )
  expect_equal(
    which(model$observations == "*"),
    c(41L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 49L, 50L, 60L, 61L, 
      62L, 63L, 64L, 65L)
  )
})
