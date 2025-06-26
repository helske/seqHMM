data("hmm_biofam")
data("mhmm_biofam")
y <- hmm_biofam$channel_names
time <- "age"
id <- "individual"
d <- stslist_to_data(
  hmm_biofam$observations, id, time, y
)
test_that("'hidden_paths' works for 'hmm'", {
  expect_error(
    out <- hidden_paths(hmm_biofam, as_stslist = TRUE),
    NA
  )
  expect_identical(
    c(table(droplevels(unlist(out)))), 
    c(`State 1` = 16075L, `State 2` = 5888L, `State 3` = 3453L, 
      `State 4` = 5094L, `State 5` = 1490L)
  )
})
test_that("'hidden_paths' works for 'mhmm'", {
  expect_error(
    out <- hidden_paths(mhmm_biofam),
    NA
  )
  expect_identical(
    c(table(out$state)), 
    c(`State 1` = 16058L, `State 2` = 5888L, `State 3` = 3438L, 
      `State 4` = 6291L, `State 5` = 238L, `State 6` = 87L)
  )
})
test_that("'hidden_paths' works for 'nhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      emission_formula = c(Marriage, Parenthood, Residence) ~ 1, n_states = 5, 
      data = d, id = id, time = time,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 1, lambda = 1, method = "DNM", check_rank = FALSE
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(
    c(table(out$state)), 
    c(`State 1` = 16075L, `State 2` = 5888L, `State 3` = 3453L, 
    `State 4` = 5094L, `State 5` = 1490L)
  )
  
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      emission_formula = Marriage ~ 1, n_states = 3,
      data = d, id = id, time = time,
      restarts = 2, maxeval = 1, method = "DNM",
      control_restart = list(maxeval = 1), check_rank = FALSE
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(
    c(table(out$state)), 
    c(`State 1` = 6692L, `State 2` = 16409L, `State 3` = 8899L)
  )
})

test_that("'hidden_paths' works for 'mnhmm'", {
  set.seed(1)
  y <- c(Marriage, Parenthood, Residence) ~ 1
  expect_error(
    fit <- estimate_mnhmm(
      emission_formula = y, n_states = 3,
      data = d, id = id, time = time, n_clusters = 2, maxeval = 1,
      maxeval_em_dnm = 1, check_rank = FALSE
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(
    c(table(out$state)), 
    c(`State 1` = 14575L, `State 2` = 1515L, `State 3` = 15910L)
  )
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      Marriage ~ 1, n_states = 3,
      data = d, id = id, time = time, n_clusters = 2,
      restarts = 2, maxeval = 1, method = "DNM",
      control_restart = list(maxeval = 1), check_rank = FALSE
    ),
    NA
  )
  expect_error(
    out <- hidden_paths(fit),
    NA
  )
  expect_identical(
    c(table(out$state)), 
    c(`State 1` = 363L, `State 2` = 9520L, `State 3` = 22117L)
  )
})
