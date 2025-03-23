
test_that("'posterior_probs' works for 'hmm'", {
  data("hmm_biofam")
  expect_error(
    out <- posterior_probs(hmm_biofam),
    NA
  )
  expect_gte(min(out$probability), 0 - 1e-12)
  expect_lte(max(out$probability), 1 + 1e-12)
})
test_that("'posterior_probs' works for 'mhmm'", {
  data("mhmm_biofam")
  expect_error(
    out <- posterior_probs(mhmm_biofam),
    NA
  )
  expect_gte(min(out$probability), 0 - 1e-12)
  expect_lte(max(out$probability), 1 + 1e-12)
})
test_that("'posterior_probs' works for 'nhmm'", {
  data("hmm_biofam")
  d <- stslist_to_data(
    hmm_biofam$observations, "id", "time", 
    hmm_biofam$channel_names
  )
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$channel_names, 
      data = d, time = "time", id = "id", n_states = 5,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 2, method = "DNM"
    ),
    NA
  )
  expect_error(
    out <- posterior_probs(fit),
    NA
  )
  expect_gte(min(out$probability), 0 - 1e-12)
  expect_lte(max(out$probability), 1 + 1e-12)
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      "Marriage", 
      data = d, time = "time", id = "id", n_states = 3,
      maxeval = 1, method = "DNM"
    ),
    NA
  )
  expect_error(
    out <- posterior_probs(fit),
    NA
  )
  expect_gte(min(out$probability), 0 - 1e-12)
  expect_lte(max(out$probability), 1 + 1e-12)
})

test_that("'posterior_probs' works for 'mnhmm'", {
  data("hmm_biofam")
  d <- stslist_to_data(
    hmm_biofam$observations, "id", "time", 
    hmm_biofam$channel_names
  )
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$channel_names, 
      data = d, time = "time", id = "id", n_states = 5, n_clusters = 2,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 2, method = "EM-DNM"
    ),
    NA
  )
  expect_error(
    out <- posterior_probs(fit),
    NA
  )
  expect_gte(min(out$probability), 0 - 1e-12)
  expect_lte(max(out$probability), 1 + 1e-12)
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      "Marriage", 
      data = d, time = "time", id = "id", n_states = 3, n_clusters = 2,
      maxeval = 1, method = "EM"
    ),
    NA
  )
  expect_error(
    out <- posterior_probs(fit),
    NA
  )
  expect_gte(min(out$probability), 0 - 1e-12)
  expect_lte(max(out$probability), 1 + 1e-12)
})
