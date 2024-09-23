test_that("'iv_X' and 'tv_X' works", {
  x1 <- array(
    c(rbind(1:4, 1:4, 1:4), rbind(1:4, 1:4, 1:4)), c(3, 4, 2)
  )
  expect_true(tv_X(x1))
  expect_false(iv_X(x1))
  x2 <- array(
    c(rbind(1:4, 1:4, 1:4), -rbind(1:4, 1:4, 1:4)), c(3, 4, 2)
  )
  expect_true(tv_X(x2))
  expect_true(iv_X(x2))
  x3 <- array(
    c(cbind(1:4, 1:4, 1:4), -cbind(1:4, 1:4, 1:4)), c(4, 3, 2)
  )
  expect_false(tv_X(x3))
  expect_true(iv_X(x3))
  x4 <- array(
    c(cbind(1:4, 1:4, 1:4), cbind(1:4, 1:4, 1:4)), c(4, 3, 2)
  )
  expect_false(tv_X(x4))
  expect_false(iv_X(x4))
})
