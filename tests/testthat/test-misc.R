test_that("check for list of formulas work ", {
  expect_true(is_formula(~1, 1L))
  expect_true(is_formula(list(~x), 1L))
  expect_true(is_formula(~x + y, 1L))
  expect_true(is_formula(list(~x, ~y), 2L))             
  expect_false(is_formula(list(~x, "not a formula")))
  expect_false(is_formula("~x + y"))
})
test_that("'cluster_names' work for 'mhmm' objects", {
  expect_equal(cluster_names(mhmm_biofam), factor(paste0("Cluster ", 1:3)))
  expect_error(cluster_names(mhmm_biofam) <- 1:3, NA)
  expect_equal(cluster_names(mhmm_biofam), factor(1:3))
})

test_that("'cluster_names' work for 'mnhmm' objects", {
  data("leaves")
  d <- leaves[leaves$workplace %in% seq_len(10), ]
  fit <- estimate_mnhmm(
    n_states = 3,
    n_clusters = 2,
    initial_formula = ~ year,
    transition_formula = ~ 1,
    emission_formula = leave ~ reform2013 + occupation,
    cluster_formula = ~ 1,
    data = d,
    id = "workplace",
    time = "father",
    method = "DNM",
    lambda = 0.1,
    maxeval = -1
  )
  expect_equal(cluster_names(fit), factor(paste0("Cluster ", 1:2)))
  expect_error(cluster_names(fit) <- 1:2, NA)
  expect_equal(cluster_names(fit), factor(1:2))
})


test_that("'state_names' work for 'hmm' objects", {
  expect_equal(
    state_names(hmm_biofam), factor(paste0("State ", 1:5))
  )
  expect_error(
    state_names(hmm_biofam) <- letters[1:5],
    NA
  )
})

test_that("'state_names' work for 'mhmm' objects", {
  expect_equal(
    state_names(mhmm_biofam), 
    list(
      `Cluster 1` = factor(paste("State", 1:4)), 
      `Cluster 2` = factor(paste("State", 1:4)), 
      `Cluster 3` = factor(paste("State", 1:6))
    )
  )
  expect_error(state_names(mhmm_biofam) <- list(1:4, 1:4, 1:6), NA)
  expect_equal(
    state_names(mhmm_biofam), 
    list(`Cluster 1` = factor(1:4), `Cluster 2` = factor(1:4), 
         `Cluster 3` = factor(1:6))
  )
})

test_that("'state_names' work for 'nhmm' objects", {
  expect_equal(
    state_names(fanhmm_leaves), factor(paste0("State ", 1:3))
  )
  expect_error(
    state_names(fanhmm_leaves) <- letters[1:3],
    NA
  )
  expect_error(
    state_names(fanhmm_leaves) <- 1,
    "Number of state names does not match with the number of states."
  )
})

test_that("'isColor works", {
  expect_equal(isColor("red"), c(red = TRUE))
  expect_equal(isColor(1), TRUE)
  expect_equal(isColor("a"), c(a = FALSE))
})

test_that("'mc_to_sc' works", {
  expect_error(
    mc_to_sc(hmm_biofam),
    NA
  )
  expect_error(
    mc_to_sc(mhmm_biofam),
    NA
  )
  expect_error(
    mc_to_sc_data(mhmm_biofam$observations),
    NA
  )
})

test_that("'most_probable_cluster' works", {
  
  expect_error(
    most_probable_cluster(mhmm_biofam, type = "abs"),
    "Argument `type` must be either \"viterbi\" or \"posterior\"."
  )
  expect_error(
    most_probable_cluster(mhmm_biofam, hp = 1),
    "Argument `hp` must be a <data.table> object from `hidden_paths\\(\\)`."
  )
  expect_equal(
    most_probable_cluster(mhmm_biofam)[7:11],
    factor(paste0("Cluster ", c(3, 2, 2, 2, 1)))
  )
  expect_equal(
    most_probable_cluster(mhmm_biofam, type = "posterior")[7:11],
    factor(paste0("Cluster ", c(3, 2, 2, 2, 1)))
  )
})
