
test_that("'forward_bacward' works for 'mhmm'", {
  data("mhmm_biofam")
  expect_error(
    fb <- forward_backward(mhmm_biofam),
    NA
  )
  expect_lte(max(fb$forward_probs), 1)
  expect_gte(min(fb$forward_probs), 0)
  expect_lte(max(fb$backward_probs), 1)
  expect_gte(min(fb$backward_probs), 0)
  expect_gte(min(fb$scaling_factors), 1)
  expect_equal(apply(fb$forward_probs,2:3,sum)

  expect_error(
    fb_log <- forward_backward(mhmm_biofam, log_space = TRUE),
    NA
  )
  expect_lte(max(fb$forward_probs), 1)
  expect_gte(min(fb$forward_probs), 0)
  expect_lte(max(fb$backward_probs), 1)
  expect_gte(min(fb$backward_probs), 0)
  expect_gte(min(fb$scaling_factors), 1)
})

