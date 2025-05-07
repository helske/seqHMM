# create test data
set.seed(123)
s <- 4
k <- 3
obs <- suppressMessages(
  seqdef(matrix(sample(letters[1:s], 50, replace = TRUE), ncol = 10))
)
obs2 <- suppressMessages(
  seqdef(matrix(sample(LETTERS[1:s], 50, replace = TRUE), ncol = 10))
)
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
    model <- build_lcm(list(y1 = obs, y2 = obs2), n_clusters = k),
    NA
  )
  expect_error(
    build_lcm(obs, n_clusters = k,
              emission_probs = cbind(1, matrix(0, 2, s - 1))),
    NA
  )
  expect_error(
    build_lcm(list(obs, obs), n_clusters = k),
    NA
  )
  expect_warning(
    model <- build_lcm(
      list(obs, obs), n_clusters = k,
      channel_names = 1:2, 
      cluster_names = letters[1:(k + 1)]),
    "The length of `cluster_names` does not match the number of clusters. Names were not used."
  )
  expect_equal(
    cluster_names(model),
    factor(paste("Class", seq_len(k)))
  )
})
test_that("build_lcm errors with incorrect dims", {
  expect_error(
    build_lcm(obs, emission_probs = diag(2)),
    "Number of columns in `emission_probs` is not equal to the number of symbols."
  )
})

test_that("build_lcm formula works", {
  expect_error(
    build_lcm(obs, n_clusters = k, formula = ~ 1),
    "If `formula` is provided, `data` must be a <data.frame> object\\."
  )
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
    factor(c("a", "b", "c", "d"))
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
