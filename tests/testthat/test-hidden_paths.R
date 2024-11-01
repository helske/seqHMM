
test_that("'hidden_paths' works for 'hmm'", {
  data("hmm_biofam")
  expect_error(
    out <- hidden_paths(hmm_biofam),
    NA
  )
  expect_identical(sum(table(unlist(out))), 32000L)
  expect_identical(length(table(unlist(out))), 7L)
})
test_that("'hidden_paths' works for 'mhmm'", {
  data("mhmm_biofam")
  expect_error(
    out <- hidden_paths(mhmm_biofam),
    NA
  )
  expect_identical(sum(table(unlist(out))), 32000L)
  expect_identical(length(table(unlist(out))), 16L)
})
test_that("'hidden_paths' works for 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations, n_states = 5,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 1
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(sum(table(unlist(out))), 32000L)
  expect_identical(length(table(unlist(out))), 7L)
  
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations[[1]], n_states = 3,
      restarts = 2, maxeval = 1
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(sum(table(unlist(out))), 32000L)
  expect_identical(length(table(unlist(out))), 5L)
})

test_that("'hidden_paths' works for 'mnhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations, n_states = 3, n_clusters = 2, maxeval = 1
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(sum(table(unlist(out))), 32000L)
  expect_identical(length(table(unlist(out))), 8L)
  
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations[[1]], n_states = 3, n_clusters = 2,
      restarts = 2, maxeval = 1
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(sum(table(unlist(out))), 32000L)
  expect_identical(length(table(unlist(out))), 8L)
})
