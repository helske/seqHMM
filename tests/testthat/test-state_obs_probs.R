
test_that("'state_obs_probs' works for multichannel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations, n_states = 5,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 1, method = "DNM"
    ),
    NA
  )
  obs <- create_obsArray(fit)
  expect_error(
    out <- state_obs_probs_nhmm_multichannel(
      fit$etas$pi, fit$X_pi, fit$etas$A, fit$X_A, 
      fit$etas$B, fit$X_B, obs, fit$sequence_lengths, 
      attr(fit$X_pi, "icpt_only"), attr(fit$X_A, "icpt_only"), 
      attr(fit$X_B, "icpt_only"), attr(fit$X_A, "iv"), 
      attr(fit$X_B, "iv"), attr(fit$X_A, "tv"), attr(fit$X_B, "tv"),
      start = 3L),
    NA
  )
  expect_gte(min(unlist(out$obs_prob)), 0)
  expect_lte(max(unlist(out$obs_prob)), 1)
  expect_gte(min(out$state_prob), 0)
  expect_lte(max(out$state_prob), 1)
  expect_true(all(abs(apply(out$obs_prob[[1]], 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$obs_prob[[2]], 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$obs_prob[[3]], 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$state_prob, 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  obs[, 3:16, ] <- fit$n_symbols
  f <- forward_nhmm_multichannel(
    fit$etas$pi, fit$X_pi,
    fit$etas$A, fit$X_A,
    fit$etas$B, fit$X_B,
    obs,
    fit$sequence_lengths, attr(fit$X_pi, "icpt_only"), 
    attr(fit$X_A, "icpt_only"), attr(fit$X_B, "icpt_only"),
    attr(fit$X_A, "iv"), attr(fit$X_B, "iv"), attr(fit$X_A, "tv"), 
    attr(fit$X_B, "tv"))
  a <- f[, , 1]
  expect_equal(
    exp(apply(a, 2, function(x) x - seqHMM:::logSumExp(x))),
    out$state_prob[, , 1]
  )
})

test_that("'state_obs_probs' works for single-channel 'nhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      hmm_biofam$observations[[1]][1:100,], n_states = 3,
      restarts = 2, maxeval = 2, lambda = 1, method = "DNM"
    ),
    NA
  )
  obs <- create_obsArray(fit)[1L, , ]
  expect_error(
    out <- state_obs_probs_nhmm_singlechannel( 
      fit$etas$pi, fit$X_pi, fit$etas$A, fit$X_A, 
      fit$etas$B, fit$X_B, obs, fit$sequence_lengths, 
      attr(fit$X_pi, "icpt_only"), attr(fit$X_A, "icpt_only"), 
      attr(fit$X_B, "icpt_only"), attr(fit$X_A, "iv"), 
      attr(fit$X_B, "iv"), attr(fit$X_A, "tv"), attr(fit$X_B, "tv"),
      start = 1L),
    NA
  )
  expect_gte(min(out$obs_prob), 0)
  expect_lte(max(out$obs_prob), 1)
  expect_gte(min(out$state_prob), 0)
  expect_lte(max(out$state_prob), 1)
  expect_true(all(abs(apply(out$obs_prob, 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$state_prob, 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
})

test_that("'state_obs_probs' works for multichannel 'mnhmm'", {
  data("hmm_biofam")
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations, n_states = 3, n_clusters = 2,
      maxeval = 1, maxeval_em_dnm = 1
    ),
    NA
  )
  obs <- create_obsArray(fit)
  expect_error(
    out <- state_obs_probs_mnhmm_multichannel( 
      fit$etas$omega, fit$X_omega, fit$etas$pi, fit$X_pi, fit$etas$A, fit$X_A, 
      unlist(fit$etas$B, recursive = FALSE), fit$X_B, obs, fit$sequence_lengths, 
      attr(fit$X_omega, "icpt_only"), attr(fit$X_pi, "icpt_only"), 
      attr(fit$X_A, "icpt_only"), attr(fit$X_B, "icpt_only"), 
      attr(fit$X_A, "iv"), attr(fit$X_B, "iv"), attr(fit$X_A, "tv"), 
      attr(fit$X_B, "tv"), start = 3L),
    NA
  )
  expect_gte(min(unlist(out$obs_prob)), 0)
  expect_lte(max(unlist(out$obs_prob)), 1)
  expect_gte(min(out$state_prob), 0)
  expect_lte(max(out$state_prob), 1)
  expect_true(all(abs(apply(out$obs_prob[[1]], 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$obs_prob[[2]], 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$obs_prob[[3]], 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$state_prob, 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
})

test_that("'state_obs_probs' works for single-channel 'mnhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      hmm_biofam$observations[[1]], n_states = 4, n_clusters = 2,
      restarts = 2, maxeval = 1, method = "DNM", algorithm = "NLOPT_LN_COBYLA"
    ),
    NA
  )
  obs <- create_obsArray(fit)[1L, , ]
  expect_error(
    out <- state_obs_probs_mnhmm_singlechannel( 
      fit$etas$omega, fit$X_omega, fit$etas$pi, fit$X_pi, fit$etas$A, fit$X_A, 
      fit$etas$B, fit$X_B, obs, fit$sequence_lengths, 
      attr(fit$X_omega, "icpt_only"), attr(fit$X_pi, "icpt_only"), 
      attr(fit$X_A, "icpt_only"), attr(fit$X_B, "icpt_only"), 
      attr(fit$X_A, "iv"), attr(fit$X_B, "iv"), attr(fit$X_A, "tv"), 
      attr(fit$X_B, "tv"), start = 3L),
    NA
  )
  expect_gte(min(out$obs_prob), 0)
  expect_lte(max(out$obs_prob), 1)
  expect_gte(min(out$state_prob), 0)
  expect_lte(max(out$state_prob), 1)
  expect_true(all(abs(apply(out$obs_prob, 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
  expect_true(all(abs(apply(out$state_prob, 2:3, sum) - 1) < sqrt(.Machine$double.eps)))
})
