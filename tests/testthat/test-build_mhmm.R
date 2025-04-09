# create test data
set.seed(123)
s <- 2:3
k <- length(s)
obs <- suppressMessages(
  seqdef(matrix(sample(letters[1:5], 50, replace = TRUE), ncol = 10))
)

test_that("build_mhmm returns object of class 'mhmm'", {
  expect_error(
    model <- build_mhmm(obs, n_states = s),
    NA
  )
  expect_s3_class(
    model,
    "mhmm"
  )
  expect_error(
    build_mhmm(obs, initial_probs = simulate_initial_probs(s, k),
              emission_probs = simulate_emission_probs(s, 5, k),
              transition_probs = simulate_transition_probs(s, k)),
    NA
  )
})
test_that("build_mhmm errors with incorrect dims", {
  expect_error(
    build_mhmm(obs, initial_probs = simulate_initial_probs(s, k+1),
               emission_probs = simulate_emission_probs(s, 5, k),
               transition_probs = simulate_transition_probs(s, k)),
    "Length of a `initial_probs` does not match with the length of `transition_probs`."
  )
  expect_error(
    build_mhmm(obs, initial_probs = simulate_initial_probs(rev(s), k),
               emission_probs = simulate_emission_probs(s, 5, k),
               transition_probs = simulate_transition_probs(s, k)),
    "Length of `initial_probs` is not equal to the number of hidden states."
  )
  expect_error(
    build_mhmm(obs, initial_probs = simulate_initial_probs(s, k),
               emission_probs = simulate_emission_probs(rev(s), 5, k),
               transition_probs = simulate_transition_probs(s, k)),
    "`emission_probs` do not sum to one."
  )
  expect_error(
    build_mhmm(obs, initial_probs = simulate_initial_probs(s, k),
               emission_probs = simulate_emission_probs(s, 4, k),
               transition_probs = simulate_transition_probs(s, k)),
    "Number of columns in `emission_probs` is not equal to the number of symbols."
  )
  expect_error(
    build_mhmm(obs, initial_probs = simulate_initial_probs(s, k),
               emission_probs = simulate_emission_probs(s, 5, k)),
    "Provide either `n_states` or all three of `initial_probs`, `transition_probs`, and `emission_probs`."
  )
})

test_that("build_mhmm formula works", {
  expect_error(
    build_mhmm(obs, n_states = k, formula = ~ 1),
    "`n_states` is of length 1, leading to an ordinary HMM. Please use `build_hmm\\(\\)` instead\\."
  )
  expect_error(
    build_mhmm(obs, n_states = s, formula = ~ 1),
    "If `formula` is provided, `data` must be a <data\\.frame> object\\."
  )
  # default error message from model.matrix, could be more informative..
  expect_error(
    build_mhmm(obs, n_states = s, formula = ~ x, data = data.frame(y = 1)),
    "object 'x' not found"
  )
  expect_error(
    build_mhmm(obs, n_states = s, formula = ~ x, data = data.frame(x = 1)),
    paste0(
      "Number of subjects in data for covariates does not match the ",
      "number of subjects in the sequence data."
    )
  )
  expect_error(
    model <- build_mhmm(obs, n_states = s, formula = ~ x,
              data = data.frame(x = 1:5),
              cluster_names = c("A", "B")),
    NA
  )
  expect_equal(
    names(model$emission_probs),
    c("A", "B")
  )
  expect_equal(
    lengths(model$emission_probs),
    c(A = 10, B = 15)
  )
  expect_equal(
    model$symbol_names,
    factor(c("a", "b", "c", "d", "e"))
  )
})

test_that("build_mhmm returns the correct number of states", {
  expect_error(
    model <- build_mhmm(obs, n_states = s),
    NA
  )
  expect_equal(
    lengths(model$initial_probs),
    c("Cluster 1" = 2, "Cluster 2" = 3)
  )
  expect_equal(
    lapply(model$emission_probs, dim),
    list("Cluster 1" = c(2, 5), "Cluster 2" = c(3, 5))
  )
})
