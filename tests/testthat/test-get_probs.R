test_that("'get_probs' works for multichannel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations, n_states = 5,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], verbose = FALSE
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit, nsim = 10),
    NA
  )
})
test_that("'get_probs' works for single-channel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations[[1]], n_states = 3,
      verbose = FALSE, iter = 1
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
})

test_that("'get_probs' works for multichannel 'mnhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations, n_states = 3, n_clusters = 2,
      verbose = FALSE, iter = 1
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
})

test_that("'get_probs' works for single-channel 'mnhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations[[1]], n_states = 4, n_clusters = 2,
      verbose = FALSE, iter = 1
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
})
