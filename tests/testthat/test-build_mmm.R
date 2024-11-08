# create test data
set.seed(123)
k <- 2
s <- 5
obs <- suppressMessages(
  seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))
)

test_that("build_mmm returns object of class 'mhmm'", {
  expect_error(
    model <- build_mmm(obs, n_clusters = k),
    NA
  )
  expect_error(
    model <- build_mmm(
      obs, 
      transition_probs = list(diag(s), diag(s)), 
      initial_probs = list(rep(1/s, s), rep(1/s, s))
    ),
    NA
  )
  expect_s3_class(
    model,
    "mhmm"
  )
  
})

test_that("build_mmm errors when neither clusters or probs are given", {
  expect_error(
    build_mmm(obs), 
    "Provide either `n_clusters` or both `initial_probs` and `transition_probs`."
  )
})
test_that("build_mmm errors with incorrect argument types", {
  expect_error(
    build_mmm(obs, transition = 1, initial = "a"), 
    "`transition_probs` is not a <list>."
  )
  expect_error(
    build_mmm(obs, transition = list(diag(2), diag(2)), initial = "a"), 
    "`initial_probs` is not a <list>."
  )
})
test_that("build_mmm errors with incorrect observations", {
  expect_error(
    build_mmm(1, n_clusters = k),
    paste0(
      "`observations` should be a <stslist> object created with ", 
      "`seqdef\\(\\)`, or a <list> of <stslist> objects in a multichannel ", 
      "case\\."
    )
  )
  expect_error(
    build_mmm(list(a = obs, b = obs), n_clusters = k),
    paste0("`build\\_mmm\\(\\)` can only be used for single-channel ",
           "sequence data \\(a <stslist> object\\)\\. Use ",
           "`mc\\_to\\_sc\\_data\\(\\)` to convert data into single-channel ", 
           "state sequences\\."
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
    stats::setNames(rep(s, k), paste("Cluster", 1:k))
  )
})
