# create test data
set.seed(123)
s <- 4
k <- 3
obs <- seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))

test_that("build_lcm returns object of class 'mhmm'", {
  expect_error(
    model <- build_lcm(obs, n_clusters = k),
    NA
  )
  expect_s3_class(
    model,
    "mhmm"
  )
  expect_error(
    build_lcm(obs, n_clusters = k,
              emission_probs = cbind(1, matrix(0, 2, s - 1))),
    NA
  )
})
test_that("build_lcm errors with incorrect dims", {
  expect_error(
    build_lcm(obs, emission_probs = diag(2)),
    "Number of columns in 'emission_probs' is not equal to the number of symbols."
  )
})

test_that("build_lcm formula works", {
  expect_error(
    build_lcm(obs, n_clusters = k, formula = ~ 1),
    "Argument 'data' is missing, but 'formula' was provided."
  )
  # default error message from model.matrix, could be more informative..
  expect_error(
    build_lcm(obs, n_clusters = k, formula = ~ x, data = data.frame(y = 1)),
    "object 'x' not found"
  )
  expect_error(
    build_lcm(obs, n_clusters = k, formula = ~ x, data = data.frame(x = 1)),
    paste0(
      "Number of subjects in data for covariates does not match the ",
      "number of subjects in the sequence data."
    )
  )
  expect_error(
    model <- build_lcm(obs, n_clusters = k, formula = ~ x,
              data = data.frame(x = 1:5),
              cluster_names = c("A", "B", "C")),
    NA
  )
  expect_equal(
    names(model$emission_probs),
    c("A", "B", "C")
  )
  expect_equal(
    lengths(model$emission_probs),
    c(A = 4, B = 4, C = 4)
  )
  expect_equal(
    model$symbol_names,
    c("a", "b", "c", "d")
  )
})

test_that("build_lcm returns the correct number of states", {
  expect_error(
    model <- build_lcm(obs, n_clusters = 2L),
    NA
  )
  expect_equal(
    unlist(model$initial_probs),
    c("Class 1.Class 1" = 1, "Class 2.Class 2" = 1)
  )
  expect_equal(
    lapply(model$emission_probs, dim),
    list("Class 1" = c(1, 4), "Class 2" = c(1, 4))
  )
})
