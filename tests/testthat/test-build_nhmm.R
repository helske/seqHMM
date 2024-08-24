# create test data
set.seed(123)
s <- 4
n_id <- 10
n_time <- 15
obs <- seqdef(matrix(sample(letters[1:s], n_id, replace = TRUE), ncol = n_time))

test_that("build_nhmm returns object of class 'nhmm'", {
  expect_error(
    model <- build_nhmm(obs, n_states = s),
    NA
  )
  expect_s3_class(
    model,
    "nhmm"
  )
  expect_error(
    build_nhmm(obs, s, initial_formula = ~ x, transition_formula = ~z, 
               emission_formula = ~ z, data0 = data.frame(x = rnorm(n_id)), 
               data = data.frame(z = rnorm(n_id * n_time))),
    NA
  )
})
test_that("build_hmm errors with missing 'n_states' argument", {
  expect_error(
    build_nhmm(obs),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})

test_that("build_nhmm errors with incorrect formulas", {
  expect_error(
    build_nhmm(obs, n_states = 3, initial_formula = 5),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_nhmm(obs, n_states = 3, transition_formula = 5),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_nhmm(obs, n_states = 3, emission_formula = 5),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})
test_that("build_hmm errors with missing data arguments", {
  expect_error(
    build_nhmm(obs, n_states = 3, initial_formula = ~ x),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_nhmm(obs, 3, initial_formula = ~ z, transition_formula = ~x,
               data0 = data.frame(rnorm(n_id))),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_nhmm(obs, 3, emission_formula = ~x),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})
test_that("build_nhmm errors with incorrect observations", {
  expect_error(
    build_nhmm(1, 2),
    paste0("Argument 'observations' should a 'stslist' object created with ",
           "'seqdef' function, or a list of such objects in case of multichannel data."
    )
  )
})
