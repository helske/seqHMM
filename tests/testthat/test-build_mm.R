# create test data
set.seed(123)
s <- 4
obs <- seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))

test_that("build_mm returns object of class 'hmm'", {
  expect_error(
    model <- build_mm(obs),
    NA
  )
  expect_s3_class(
    model,
    "hmm"
  )
})
test_that("build_mm errors with incorrect observations", {
  expect_error(
    build_mm(1),
    paste0(
      "Argument 'observations' should a 'stslist' object created with ",
      "'seqdef' function, or a list of such objects in case of multichannel ",
      "data."
    )
  )
  expect_error(
    build_mm(list(a = obs, b = obs)),
    paste0("The 'build_mm' function can only be used for single-channel ",
           "sequence data \\(as an stslist object\\). Use the 'mc_to_sc_data' function ",
           "to convert data into single-channel state sequences."
    )
  )
})

test_that("build_mm returns the correct number of states", {
  expect_error(
    model <- build_mm(obs),
    NA
  )
  expect_equal(
    length(model$initial_probs),
    s
  )
  expect_equal(
    dim(model$transition_probs),
    c(s, s)
  )
  expect_equal(
    dim(model$emission_probs),
    c(s, s)
  )
})

test_that("build_mm returns the correct probabilities", {
  model <- build_mm(obs)
  expect_equal(
    model$initial_probs,
    c(prop.table(table(obs[, 1]))[1:s])
  )
  expect_equal(
    rowSums(model$transition_probs),
    setNames(rep(1, s), letters[1:s])
  )
  expect_equal(
    rowSums(model$emission_probs),
    setNames(rep(1, s), letters[1:s])
  )
  expect_true(all(model$transition_probs >= 0))
  expect_true(all(model$emission_probs >= 0))
})
