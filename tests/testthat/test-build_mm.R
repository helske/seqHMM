# create test data
set.seed(123)
s <- 4
obs <- suppressMessages(
  seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))
)

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
      "`observations` should be a <stslist> object created with ", 
      "`seqdef\\(\\)`, or a <list> of <stslist> objects in a multichannel ", 
      "case\\."
    )
  )
  expect_error(
    build_mm(list(a = obs, b = obs)),
    paste0("`build\\_mm\\(\\)` can only be used for single-channel ",
           "sequence data \\(a <stslist> object\\)\\. Use ",
           "`mc\\_to\\_sc\\_data\\(\\)` to convert data into single-channel state ", 
           "sequences\\."
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
