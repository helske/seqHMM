# create test data
set.seed(123)
s <- 4
d <- 3
n_id <- 10
n_time <- 15
obs <- seqdef(matrix(sample(letters[1:s], n_id, replace = TRUE), ncol = n_time))

test_that("build_mnhmm returns object of class 'mnhmm'", {
  expect_error(
    model <- build_mnhmm(obs, n_states = s, n_clusters = d),
    NA
  )
  expect_s3_class(
    model,
    "nhmm"
  )
  expect_error(
    build_mnhmm(obs, s, d, initial_formula = ~ x, transition_formula = ~z,
                emission_formula = ~ z, cluster_formula = ~ x,
                data0 = data.frame(x = rnorm(n_id), 
                                   y = factor(sample(1:3, n_id, TRUE))), 
                data = data.frame(z = rnorm(n_id * n_time))),
    NA
  )
})
test_that("build_mhmm errors with missing 'n_states' argument", {
  expect_error(
    build_mnhmm(obs),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})
test_that("build_mhmm errors with missing 'n_clusters' argument", {
  expect_error(
    build_mnhmm(obs, s),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})
test_that("build_mnhmm errors with incorrect formulas", {
  expect_error(
    build_mnhmm(obs, n_states = 3, initial_formula = 5),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_mnhmm(obs, n_states = 3, transition_formula = 5),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_mnhmm(obs, n_states = 3, emission_formula = 5),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_mnhmm(obs, n_states = 3, cluster_formula = 5),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})
test_that("build_mnhmm errors with missing data arguments", {
  expect_error(
    build_mnhmm(obs, s, d, initial_formula = ~ x),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_mnhmm(obs, 3, 2, initial_formula = ~ z, transition_formula = ~x,
               data0 = data.frame(rnorm(n_id))),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
  expect_error(
    build_mnhmm(obs, 3, 5, emission_formula = ~x),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})
test_that("build_hmm errors with incorrect observations", {
  expect_error(
    build_mnhmm(list(1,"a"), 2),
    paste0("Argument 'observations' should a 'stslist' object created with ",
           "'seqdef' function, or a list of such objects in case of multichannel data."
    )
  )
})
