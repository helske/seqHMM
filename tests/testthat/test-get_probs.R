test_that("'get_probs' and 'coef' works for multichannel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations, n_states = 5,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 1, method = "DNM"
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
  expect_error(
    coef(fit),
    NA
  )
})
test_that("'get_probs' and 'coef' works for single-channel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations[[1]], n_states = 3, maxeval = 1, method = "DNM"
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
  expect_error(
    coef(fit),
    NA
  )
})

test_that("'get_probs' and 'coef' works for multichannel 'mnhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations, n_states = 3, n_clusters = 2, maxeval = 1,
      em_dnm_maxeval = 1
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
  expect_error(
    coef(fit),
    NA
  )
})

test_that("'get_probs' and 'coef' works for single-channel 'mnhmm'", {
  set.seed(1)
  d <- data.frame(
    group = rep(1:50, each = 16),
    time = 1:16,
    z = rnorm(16 * 50),
    w = 1:16
  )
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations[[1]][1:50, ], n_states = 4, n_clusters = 2, 
      initial_formula = ~ z, cluster_formula = ~ z, 
      transition_formula = ~w, emission_formula = ~ w, 
      data = d, time = "time", id = "group", maxeval = 1,
      em_dnm_maxeval = 1
    ),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
  expect_error(
    coef(fit),
    NA
  )
})
