
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
      emission_formula = c(Marriage, Parenthood) ~ 1, 
      data = d, time = "time", id = "id", n_states = 5,
      maxeval = 2, method = "DNM", check_rank = FALSE
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
      Marriage ~ 1, 
      data = d, time = "time", id = "id", n_states = 3,
      maxeval = 1, method = "DNM", check_rank = FALSE
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
    fit <- estimate_mnhmm(
      emission = lapply(
        hmm_biofam$channel_names, \(y) stats::as.formula(paste0(y, " ~ 1"))
      ), 
      data = d, time = "time", id = "id", n_states = 2, n_clusters = 2,
      maxeval = 2, method = "EM-DNM", check_rank = FALSE
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
      Marriage ~ 1, 
      data = d, time = "time", id = "id", n_states = 3, n_clusters = 2,
      maxeval = 1, method = "EM", check_rank = FALSE
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
