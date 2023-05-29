
test_that("'forward_bacward' works for 'mhmm'", {
  data("mhmm_biofam")
  expect_error(
    fb <- forward_backward(mhmm_biofam),
    NA
  )
  expect_gte(min(fb$forward_probs), 0)
  expect_gte(min(fb$backward_probs), 0)

  expect_error(
    fb_log <- forward_backward(mhmm_biofam, log_space = TRUE),
    NA
  )
  expect_lte(max(exp(fb_log$forward_probs)), 1)
  expect_gte(min(exp(fb_log$forward_probs)), 0)
  expect_lte(max(exp(fb_log$backward_probs)), 1)
  expect_gte(min(exp(fb_log$backward_probs)), 0)
})

