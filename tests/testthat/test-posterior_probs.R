
test_that("'posterior_probs' works for 'hmm'", {
  data("hmm_biofam")
  expect_error(
    out <- posterior_probs(
      hmm_biofam, log_space = FALSE, as_data_frame = FALSE
    ),
    NA
  )
  expect_gte(min(out), -1e-8)
  expect_lte(max(out), 1 + 1e-8)
  expect_error(
    out <- posterior_probs(hmm_biofam, log_space = TRUE),
    NA
  )
  expect_gte(min(out$probability), -1e-8)
  expect_lte(max(out$probability), 1 + 1e-8)
})
test_that("'posterior_probs' works for 'mhmm'", {
  data("mhmm_biofam")
  expect_error(
    out <- posterior_probs(
      mhmm_biofam, log_space = FALSE, as_data_frame = FALSE
    ),
    NA
  )
  expect_gte(min(out), -1e-8)
  expect_lte(max(out), 1 + 1e-8)
  expect_error(
    out <- posterior_probs(hmm_biofam, log_space = TRUE),
    NA
  )
  expect_gte(min(out$probability), -1e-8)
  expect_lte(max(out$probability), 1 + 1e-8)
})
test_that("'posterior_probs' works for 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations, n_states = 5,
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
  expect_gte(min(out$probability), -1e-8)
  expect_lte(max(out$probability), 1 + 1e-8)
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations[[1]], n_states = 3,
      maxeval = 1, method = "DNM"
    ),
    NA
  )
  expect_error(
    out <- posterior_probs(fit),
    NA
  )
  expect_gte(min(out$probability), -1e-8)
  expect_lte(max(out$probability), 1 + 1e-8)
})

test_that("'posterior_probs' works for 'mnhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations, n_states = 3, n_clusters = 2, maxeval = 1,
      method = "EM"
    ),
    NA
  )
  expect_error(
    out <- posterior_probs(fit),
    NA
  )
  expect_gte(min(out$probability), -1e-8)
  expect_lte(max(out$probability), 1 + 1e-8)
  
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations[[1]], n_states = 3, n_clusters = 2, 
      maxeval = 2, maxeval_em_dnm = 2, lambda = 1e-2
    ),
    NA
  )
  expect_error(
    out <- posterior_probs(fit),
    NA
  )
  expect_gte(min(out$probability), -1e-8)
  expect_lte(max(out$probability), 1 + 1e-8)
})
