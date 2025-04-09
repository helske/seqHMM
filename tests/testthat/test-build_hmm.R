# create test data
set.seed(123)
s <- 4
obs <- suppressMessages(
  seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))
)

test_that("build_hmm returns object of class 'hmm'", {
  expect_error(
    model <- build_hmm(obs, n_states = s),
    NA
  )
  expect_s3_class(
    model,
    "hmm"
  )
  expect_error(
    build_hmm(obs, initial_probs = c(1, 0),
              transition_probs = diag(2),
              emission_probs = cbind(1, matrix(0, 2, s - 1))),
    NA
  )
})
test_that("build_hmm errors with incorrect dims", {
  expect_error(
    build_hmm(obs, initial_probs = c(1, 0),
              transition_probs = diag(2),
              emission_probs = diag(2)),
    "Number of columns in `emission_probs` is not equal to the number of symbols."
  )
  expect_error(
    build_hmm(obs, initial_probs = c(1, 0, 0),
              transition_probs = diag(2),
              emission_probs = cbind(1, matrix(0, 2, s - 1))),
    "Length of `initial_probs` is not equal to the number of hidden states."
  )
  expect_error(
    build_hmm(obs, initial_probs = c(1, 0, 0),
              transition_probs = diag(3),
              emission_probs = cbind(1, matrix(0, 2, s - 1))),
    "`emission_probs` do not sum to one."
  )
})

test_that("build_hmm errors with incorrect observations", {
  expect_error(
    build_hmm(1, initial_probs = c(1, 0),
              transition_probs = diag(2),
              emission_probs = diag(2)),
    paste0("`observations` should be a <stslist> object created with `seqdef\\(\\)`,",
           " or a <list> of <stslist> objects in a multichannel case\\.")
  )
})

test_that("build_hmm returns the correct number of states", {
  expect_error(
    model <- build_hmm(obs, n_states = s),
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

test_that("build_hmm returns the correct probabilities", {
  model <- build_hmm(obs, n_states = s)
  expect_true(
    all(model$initial_probs >= 0)
  )
  expect_true(
    all(model$initial_probs <= 1)
  )
  expect_equal(sum(model$initial_probs), 1)
  
  expect_equal(
    rowSums(model$transition_probs),
    stats::setNames(rep(1, s), paste("State", 1:s))
  )
  expect_equal(
    rowSums(model$emission_probs),
    stats::setNames(rep(1, s), paste("State", 1:s))
  )
  expect_equal(colnames(model$emission_probs), letters[1:s])
  expect_true(all(model$transition_probs >= 0))
  expect_true(all(model$emission_probs >= 0))
  expect_true(all(model$transition_probs <= 1))
  expect_true(all(model$emission_probs <= 1))
})
