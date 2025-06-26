data("hmm_biofam")
responses <- hmm_biofam$channel_names
time <- "age"
id <- "individual"
d <- stslist_to_data(
  hmm_biofam$observations, id, time, responses
)
test_that("'forward_backward' works for 'hmm'", {
  expect_error(
    fb <- forward_backward(hmm_biofam),
    NA
  )
  expect_lte(max(fb$log_alpha), 0)
  expect_lte(max(fb$log_beta), 0)
})
test_that("'forward_backward' works for 'mhmm'", {
  expect_error(
    fb <- forward_backward(mhmm_biofam),
    NA
  )
  expect_lte(max(fb$log_alpha), 0)
  expect_lte(max(fb$log_beta), 0)
})
test_that("'forward_backward' works for multichannel 'nhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      5, c(Marriage, Parenthood, Residence) ~ 1,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 1, method = "DNM",
      data = d, id = id, time = time, check_rank = FALSE
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit),
    NA
  )
  expect_lte(max(fb$log_alpha), 0)
  expect_lte(max(fb$log_beta), 0)
})
test_that("'forward_backward' works for single-channel 'nhmm'", {
 
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      n_states = 3, Parenthood ~ 1,
      data = d, id = id, time = time,
      restarts = 2, maxeval = 2, lambda = 0.1, method = "EM",
      control_restart = list(maxeval = 2), check_rank = FALSE
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit),
    NA
  )
  expect_lte(max(fb$log_alpha), 0)
  expect_lte(max(fb$log_beta), 0)
})

test_that("'forward_backward' works for multichannel 'mnhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      c(Marriage, Parenthood, Residence) ~ 1, n_states = 3, n_clusters = 2,
      maxeval = 1, method = "EM", data = d, id = id, time = time, 
      check_rank = FALSE
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit),
    NA
  )
  expect_lte(max(fb$log_alpha), 0)
  expect_lte(max(fb$log_beta), 0)
})

test_that("'forward_backward' works for single-channel 'mnhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      emission_formula =  Marriage ~ Residence, n_states = 2, n_clusters = 2,
      data = d, id = id, time = time, maxeval = -1, check_rank = FALSE
    ),
    NA
  )
  expect_error(
    fb <- forward_backward(fit),
    NA
  )
  expect_lte(max(fb$log_alpha), 0)
  expect_lte(max(fb$log_beta), 0)
})
