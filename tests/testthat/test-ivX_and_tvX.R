test_that("'iv_X' and 'tv_X' works", {
  id <- rep(1:2, each = 5)
  time <- rep(1:5, 2)
  
  # x varies by time but not by id
  X1 <- matrix(1:5, nrow = 10, ncol = 1)
  expect_true(tv_X(X1, 1, id))
  expect_false(iv_X(X1, 1, time))
  
  # x varies by time and id
  X2 <- matrix(1:10, nrow = 10, ncol = 1)
  expect_true(tv_X(X2, 1, id))
  expect_true(iv_X(X2, 1, time))
  
  # x varies by id but not by time
  X3 <- matrix(rep(1:2, each = 5), nrow = 10, ncol = 1)
  expect_false(tv_X(X3, 1, id))
  expect_true(iv_X(X3, 1, time))
  
  # x does not vary
  X4 <- matrix(1, nrow = 10, ncol = 1)
  expect_false(tv_X(X4, 1, id))
  expect_false(iv_X(X4, 1, time))
  
  # Multiple columns, second varies by time but not by id
  X5 <- cbind(rep(1, 10), 1:5)
  expect_true(tv_X(X5, 1:2, id))
  expect_false(iv_X(X5, 1:2, time))
})

test_that("'iv_X' and 'tv_X' works with varying sequence lengths", {
  id <- c(1, 1, 1, 2, 2)
  time <- c(1, 2, 3, 1, 2)
  
  # x varies by both
  X1 <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  expect_true(tv_X(X1, 1, id))
  expect_true(iv_X(X1, 1, time))  # different values at same time points
  
  # x varies by id
  X2 <- matrix(c(1, 1, 1, 2, 2), nrow = 5, ncol = 1)
  expect_false(tv_X(X2, 1, id))
  expect_true(iv_X(X2, 1, time))
  
  # x varies by time but not by id (for common time points)
  X3 <- matrix(c(1, 2, 3, 1, 2), nrow = 5, ncol = 1)
  expect_true(tv_X(X3, 1, id))
  expect_false(iv_X(X3, 1, time))
  
  # x constant
  X4 <- matrix(5, nrow = 5, ncol = 1)
  expect_false(tv_X(X4, 1, id))
  expect_false(iv_X(X4, 1, time))
})

test_that("model_matrix_transition_formula works with intercept only", {
  data <- data.table(
    id = c(1, 2, 2, 2, 3, 3),
    time = c(1, 1:3, 1:2),
    y = sample(c("a", "b"), 6, replace = TRUE)
  )
  
  X <- model_matrix_transition_formula(
    formula = ~ 1,
    data = data,
    sequence_lengths = c(1, 3, 2),
    id_var = "id",
    time_var = "time",
    check = FALSE
  )
  
  expect_true(is.list(X))
  expect_length(X, 3)
  expect_equal(dim(X[[1]]), c(1, 1))
  expect_equal(dim(X[[2]]), c(1, 3))
  expect_equal(dim(X[[3]]), c(1, 2))
  expect_true(all(sapply(X, \(x) all(x == 1))))
  expect_equal(attr(X, "coef_names"), "(Intercept)")
  expect_true(attr(X, "icpt_only"))
  expect_false(attr(X, "tv"))
  expect_true(attr(X, "iv"))
})

test_that("model_matrix_transition_formula works with covariates", {
  data <- data.table(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    y = sample(c("a", "b"), 12, replace = TRUE),
    x1 = rnorm(12),
    x2 = rep(c(0, 1), 6)
  )
  
  X <- model_matrix_transition_formula(
    formula = ~ x1 + x2,
    data = data,
    sequence_lengths = rep(4, 3),
    id_var = "id",
    time_var = "time",
    check = FALSE
  )
  
  expect_true(is.list(X))
  expect_length(X, 3)
  expect_equal(dim(X[[1]]), c(3, 4))
  expect_equal(dim(X[[2]]), c(3, 4))
  expect_equal(dim(X[[3]]), c(3, 4))
  expect_equal(attr(X, "coef_names"), c("(Intercept)", "x1", "x2"))
  expect_false(attr(X, "icpt_only"))
  expect_true(attr(X, "tv"))  # x1 varies within id
  expect_true(attr(X, "iv"))  # x1 varies across ids at same time
})

test_that("model_matrix_transition_formula splits correctly by id", {
  set.seed(123)
  data <- data.table(
    id = rep(1:2, each = 3),
    time = rep(1:3, 2),
    y = sample(c("a", "b"), 6, replace = TRUE),
    x = c(1, 2, 3, 4, 5, 6)
  )
  
  X <- model_matrix_transition_formula(
    formula = ~ x,
    data = data,
    sequence_lengths = rep(3, 2),
    id_var = "id",
    time_var = "time",
    scale = FALSE,
    check = FALSE
  )
  expect_equal(X[[1]][2, ], c(1, 2, 3), ignore_attr = TRUE)
  expect_equal(X[[2]][2, ], c(4, 5, 6), ignore_attr = TRUE)
  expect_equal(X[[1]][1, ], c(1, 1, 1), ignore_attr = TRUE)
  expect_equal(X[[2]][1, ], c(1, 1, 1), ignore_attr = TRUE)
})

test_that("model_matrix_transition_formula handles varying sequence lengths", {
  data <- data.table(
    id = c(1, 1, 1, 2, 2),
    time = c(1, 2, 3, 1, 2),
    y = sample(c("a", "b"), 5, replace = TRUE),
    x = 1:5
  )
  
  X <- model_matrix_transition_formula(
    formula = ~ x,
    data = data,
    sequence_lengths = c(3, 2),
    id_var = "id",
    time_var = "time",
    scale = FALSE,
    check = FALSE
  )
  
  expect_equal(dim(X[[1]]), c(2, 3))
  expect_equal(dim(X[[2]]), c(2, 2))
  expect_true(attr(X, "iv"))
  expect_true(attr(X, "tv"))
})
