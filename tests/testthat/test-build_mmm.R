# create test data
set.seed(123)
k <- 2
s <- 5
obs <- seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))

test_that("build_mmm returns object of class 'mhmm'", {
  expect_error(
    model <- build_mmm(obs, n_clusters = k),
    NA
  )
  expect_s3_class(
    model,
    "mhmm"
  )
})
test_that("build_mmm errors with incorrect observations", {
  expect_error(
    build_mmm(1, n_clusters = k),
    paste0(
      "Argument 'observations' should a 'stslist' object created with ",
      "'seqdef' function, or a list of such objects in case of multichannel ",
      "data."
    )
  )
  expect_error(
    build_mmm(list(a = obs, b = obs), n_clusters = k),
    paste0("The 'build_mmm' function can only be used for single-channel ",
           "sequence data \\(as an stslist object\\). Use the 'mc_to_sc_data' function ",
           "to convert data into single-channel state sequences."
    )
  )
})

test_that("build_mmm returns the correct number of states", {
  expect_error(
    model <- build_mmm(obs, n_clusters = k),
    NA
  )
  expect_equal(
    lengths(model$initial_probs),
    setNames(rep(s, k), paste("Cluster", 1:k))
  )
})
