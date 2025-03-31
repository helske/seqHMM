test_that("simulate_hmm throws error with incorrect arguments", {
  expect_error(
    simulate_hmm(),
   "`emission_probs` must be a matrix or a list of matrices.")
})
