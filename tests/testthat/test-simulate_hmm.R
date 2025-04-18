test_that("simulate_hmm throws error with incorrect arguments", {
  expect_error(
    simulate_hmm(),
    "`emission_probs` must be a matrix or a list of matrices."
  )
  
  emission_probs <- matrix(c(0.5, 0.2, 0.5, 0.8), 2, 2)
  transition_probs <- matrix(c(5 / 6, 1 / 6, 1 / 6, 5 / 6), 2, 2)
  initial_probs <- c(1, 0)
  set.seed(1)
  
  expect_error(
    sim <- simulate_hmm(
      n_sequences = 3, initial_probs = initial_probs,
      transition_probs = transition_probs,
      emission_probs = emission_probs,
      sequence_length = 10
    ),
    NA
  )
  expect_equal(names(sim), c("observations", "states"))
  expect_true(is_stslist(sim, 2L))
})
