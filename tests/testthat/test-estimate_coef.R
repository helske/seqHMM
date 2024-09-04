
test_that("'estimate_coef' works for 'mhmm'", {
  set.seed(123)
  data("mhmm_biofam")
  expect_error(
    out <- estimate_coef(mhmm_biofam),
    NA
  )
  expect_gte(out, logLik(mhmm_biofam))
})

test_that("'estimate_coef' errors for 'hmm'", {
  set.seed(123)
  data("hmm_biofam")
  expect_error(
    out <- estimate_coef(hmm_biofam),
    "`model` must be a `mhmm` object."
  )
})
test_that("'estimate_coef' errors for incorrect 'threads'", {
  data("mhmm_biofam")
  expect_error(
    out <- estimate_coef(mhmm_biofam, threads = "a"),
    "Argument `threads` must be a single positive integer"
  )
})
