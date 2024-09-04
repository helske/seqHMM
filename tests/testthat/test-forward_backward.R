
test_that("'forward_backward' works for 'hmm'", {
  data("hmm_biofam")
  expect_error(
    fb <- forward_backward(hmm_biofam, log_space = FALSE),
    NA
  )
  expect_gte(min(fb$forward_probs), 0)
  expect_gte(min(fb$backward_probs), 0)
  
  expect_error(
    fb_log <- forward_backward(mhmm_biofam),
    NA
  )
  expect_lte(max(exp(fb_log$log_probability)), 1)
  expect_gte(min(exp(fb_log$log_probability)), 0)
})
test_that("'forward_backward' works for 'mhmm'", {
  data("mhmm_biofam")
  expect_error(
    fb <- forward_backward(mhmm_biofam, log_space = FALSE),
    NA
  )
  expect_gte(min(fb$forward_probs), 0)
  expect_gte(min(fb$backward_probs), 0)
  
  expect_error(
    fb_log <- forward_backward(mhmm_biofam),
    NA
  )
  expect_lte(max(exp(fb_log$log_probability)), 1)
  expect_gte(min(exp(fb_log$log_probability)), 0)
})
test_that("'forward_backward' works for multichannel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations, n_states = 5,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ],
      iter = 1, verbose = FALSE, hessian = FALSE,
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit, as_data_frame = FALSE),
    NA
  )
  expect_gte(min(fb$forward_probs), -200)
  expect_gte(min(fb$backward_probs), -200)
  expect_lte(max(fb$forward_probs), 0)
  expect_lte(max(fb$backward_probs), 0)
  
  expect_error(
    fb_log <- forward_backward(mhmm_biofam),
    NA
  )
  expect_lte(max(exp(fb_log$log_probability)), 1)
  expect_gte(min(exp(fb_log$log_probability)), 0)
})
test_that("'forward_backward' works for single-channel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations[[1]], n_states = 3,
      iter = 1, verbose = FALSE, restarts = 2, threads = 1, hessian = FALSE
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit, as_data_frame = FALSE),
    NA
  )
  expect_gte(min(fb$forward_probs), -60)
  expect_gte(min(fb$backward_probs), -60)
  expect_lte(max(fb$forward_probs), 0)
  expect_lte(max(fb$backward_probs), 0)
  
  expect_error(
    fb_log <- forward_backward(fit),
    NA
  )
  expect_lte(max(exp(fb_log$log_probability)), 1)
  expect_gte(min(exp(fb_log$log_probability)), 0)
})

test_that("'forward_backward' works for multichannel 'mnhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations, n_states = 3, n_clusters = 2,
      iter = 1, verbose = FALSE, hessian = FALSE,
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit, as_data_frame = FALSE),
    NA
  )
  expect_gte(min(fb$forward_probs), -130)
  expect_gte(min(fb$backward_probs), -120)
  expect_lte(max(fb$forward_probs), 0)
  expect_lte(max(fb$backward_probs), 0)
  
  expect_error(
    fb_log <- forward_backward(fit),
    NA
  )
  expect_lte(max(exp(fb_log$log_probability)), 1)
  expect_gte(min(exp(fb_log$log_probability)), 0)
})

test_that("'forward_backward' works for single-channel 'mnhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations[[1]], n_states = 4, n_clusters = 2,
      iter = 1, verbose = FALSE, restarts = 2, threads = 1, hessian = FALSE
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit, as_data_frame = FALSE),
    NA
  )
  expect_gte(min(fb$forward_probs), -50)
  expect_gte(min(fb$backward_probs), -50)
  expect_lte(max(fb$forward_probs), 0)
  expect_lte(max(fb$backward_probs), 0)
  
  expect_error(
    fb_log <- forward_backward(fit),
    NA
  )
  expect_lte(max(exp(fb_log$log_probability)), 1)
  expect_gte(min(exp(fb_log$log_probability)), 0)
})
