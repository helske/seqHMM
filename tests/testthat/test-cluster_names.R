
test_that("'cluster_names' works for 'mhmm'", {
  expect_equal(
    cluster_names(mhmm_biofam),
    factor(c("Cluster 1", "Cluster 2", "Cluster 3"))
  )
  expect_error(
    cluster_names(mhmm_biofam) <- letters[1:3],
    NA
  )
  expect_equal(
    cluster_names(mhmm_biofam),
    factor(letters[1:3])
  )
})

test_that("'cluster_names.mhmm' errors", {
  expect_error(
    cluster_names(mhmm_biofam) <- 1:4,
    "New cluster names should be a vector of length 3."
  )
  expect_error(
    cluster_names(mhmm_biofam) <- diag(2),
    "New cluster names should be a vector of length 3."
  )
})
