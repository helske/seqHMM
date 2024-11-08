test_that("build_mm returns object of class 'hmm'", {
  set.seed(123)
  s <- 4
  obs_matrix <- matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10)
  obs_matrix[1:3, 10] <- NA
  obs_matrix[5, 5] <- NA
  obs_matrix[4, 0] <- "z"
  obs_matrix[5, 10] <- "z"
  obs <- suppressMessages(seqdef(obs_matrix))
  expect_message(
    model <- build_mm(obs),
    "Sequences contain missing values, initial and transition probabilities estimated via EM."
  )
  expect_s3_class(
    model,
    "hmm"
  )
  set.seed(123)
  s <- 4
  obs_matrix <- matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10)
  obs_matrix[4, 10] <- "z"
  obs_matrix[5, 10] <- "z"
  obs <- suppressMessages(seqdef(obs_matrix))
  expect_warning(
    model <- build_mm(obs),
    "There are no observed transitions from some of the symbols."
  )
})
set.seed(123)
s <- 4
obs_matrix <- matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10)
obs <- suppressMessages(seqdef(obs_matrix))
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
  set.seed(123)
  s <- 4 
  obs <- suppressMessages(
    seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))
  )
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
  set.seed(123)
  s <- 4 
  obs <- suppressMessages(
    seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))
  )
  model <- build_mm(obs)
  expect_equal(
    model$initial_probs,
    c(prop.table(table(obs[, 1]))[1:s])
  )
  expect_equal(
    rowSums(model$transition_probs),
    stats::setNames(rep(1, s), letters[1:s])
  )
  expect_equal(
    rowSums(model$emission_probs),
    stats::setNames(rep(1, s), letters[1:s])
  )
  expect_true(all(model$transition_probs >= 0))
  expect_true(all(model$emission_probs >= 0))
})
