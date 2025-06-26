test_that("'predict' works 'fanhmm'", {
  data("fanhmm_leaves")
  data("leaves")
  d <- leaves[leaves$workplace %in% seq_len(10), ]
  expect_error(
    out <- predict(fanhmm_leaves, newdata = d),
    NA
  )
  expect_error(
    predict(fanhmm_leaves),
    "Argument `newdata` must be a <data.frame> object."
  )
  expect_error(
    predict(fanhmm_leaves, newdata = d, newdata2 = "a"),
    "Argument `newdata2` must be a <data.frame> object."
  )
  expect_error(
    predict(fanhmm_leaves, newdata = d, type = "abc"),
    'Argument `type` must be \"state\", \"response\", \"transition\", \"emission\", or a combination of these.'
  )
})
test_that("'predict' works multichannel 'nhmm'", {
  data("hmm_biofam")
  y <- hmm_biofam$channel_names
  age <- "age"
  id_var <- "individual"
  d <- stslist_to_data(
    hmm_biofam$observations, id_var, age, y
  )
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      n_states = 3, emission_formula = c(Marriage, Parenthood) ~ Residence,
      data = d, time = age, id = id_var,
      maxeval = 3, method = "DNM", check_rank = FALSE
    ),
    NA
  )
  fit <- bootstrap_coefs(fit, nsim = 5)
  d$Marriage[d$age > 25] <- NA
  d$Parenthood[d$age > 28] <- NA
  expect_error(
    out <- predict(
      fit, newdata = d, condition = "Residence", probs = c(0.1, 0.9)
    ),
    NA
  )
})
test_that("'predict' works 'mnhmm'", {
  data("leaves")
  d <- leaves[leaves$workplace %in% seq_len(10), ]
  fit <- estimate_mnhmm(
    n_states = 2,
    n_clusters = 2,
    initial_formula = ~ year,
    transition_formula = ~ 1,
    emission_formula = leave ~ reform2013 + occupation,
    cluster_formula = ~ 1,
    data = d,
    id = "workplace",
    time = "father",
    method = "DNM",
    lambda = 0.1
  )
  expect_equal(fit$estimation_results$return_code, 3)
  fit <- bootstrap_coefs(fit, nsim = 5)
  expect_error(
    out <- predict(fit, newdata = d, probs = c(0.1, 0.9)),
    NA
  )
})
