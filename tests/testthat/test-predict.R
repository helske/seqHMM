test_that("'predict' works 'fanhmm'", {
  data("fanhmm_leaves")
  data("leaves")
  d <- leaves[leaves$workplace %in% seq_len(10), ]
  expect_error(
    out <- predict(fanhmm_leaves, newdata = d),
    NA
  )
  expect_error(
    predict(fanhmm_leaves),
    "Argument `newdata` must be a <data.frame> object."
  )
  expect_error(
    predict(fanhmm_leaves, newdata = d, newdata2 = "a"),
    "Argument `newdata2` must be a <data.frame> object."
  )
  expect_error(
    predict(fanhmm_leaves, newdata = d, type = "abc"),
    'Argument `type` must be \"observations\", \"states\", \"conditionals\", or a combination of these.'
  )
})
test_that("'predict' works 'mnhmm'", {
  data("leaves")
  d <- leaves[leaves$workplace %in% seq_len(10), ]
  fit <- estimate_mnhmm(
    n_states = 2,
    n_clusters = 2,
    initial_formula = ~ year,
    transition_formula = ~ 1,
    emission_formula = leave ~ reform2013 + occupation,
    cluster_formula = ~ 1,
    data = d,
    id = "workplace",
    time = "father",
    method = "DNM",
    lambda = 0.1
  )
  expect_equal(fit$estimation_results$return_code, 3)
  expect_error(
    out <- predict(fit, newdata = d),
    "no applicable method for 'predict' applied to an object of class \"mnhmm\""
  )
})
