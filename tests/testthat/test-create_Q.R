test_that("Test that QR decomposition of A works", {
  create_QR <- function(K) {
    A <- diag(K)
    A[K, 1:(K - 1)] <- -1
    A[K, K] <- 0
    qr.Q(qr(A))[, 1:(K - 1)]
  }
  set.seed(123)
  S <- 10
  x <- stats::rnorm(S - 1)
  expect_equal(create_QR(S), -create_Q(S))
  expect_equal(create_QR(S) %*% x, -create_Q(S) %*% x)
})

