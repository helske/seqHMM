data("hmm_biofam")
y <- hmm_biofam$channel_names
age <- "age"
id_var <- "individual"
d <- stslist_to_data(
  hmm_biofam$observations, id_var, age, y
)
test_that("'get_probs' and 'coef' works for multichannel 'nhmm'", {
  
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      y, n_states = 5,
      data = d, time = age, id = id_var,
      inits = hmm_biofam[
        c("initial_probs", "transition_probs", "emission_probs")
      ], maxeval = 1, method = "DNM"
    ),
    NA
  )
  expect_error(
    p <- get_initial_probs(fit),
    NA
  )
  expect_error(
    p <- get_emission_probs(fit),
    NA
  )
  expect_error(
    p <- get_transition_probs(fit),
    NA
  )
  expect_error(
    coef(fit),
    NA
  )
  expect_error(
    marginals <- get_marginals(fit),
    NA
  )
})
test_that("'get_probs' and 'coef' works for single-channel 'nhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_nhmm(
      y[1], n_states = 3, data = d, time = age, id = id_var, 
      maxeval = 1, method = "DNM"
    ),
    NA
  )
  expect_error(
    p <- get_initial_probs(fit),
    NA
  )
  expect_error(
    p <- get_emission_probs(fit),
    NA
  )
  expect_error(
    p <- get_transition_probs(fit),
    NA
  )
  expect_error(
    coef(fit),
    NA
  )
})

test_that("'get_probs' and 'coef' works for multichannel 'mnhmm'", {
  set.seed(1)
  expect_error(
    fit <- estimate_mnhmm(
      y, data = d, time = age, id = id_var, 
      n_states = 3, n_clusters = 2, maxeval = 1,
      maxeval_em_dnm = 1
    ),
    NA
  )
  expect_error(
    p <- get_initial_probs(fit),
    NA
  )
  expect_error(
    p <- get_emission_probs(fit),
    NA
  )
  expect_error(
    p <- get_transition_probs(fit),
    NA
  )
  expect_error(
    coef(fit),
    NA
  )
})

test_that("'get_probs' and 'coef' works for single-channel 'mnhmm'", {
  set.seed(1)
  d <- data.frame(
    group = rep(1:50, each = 16),
    time = 1:16,
    z = stats::rnorm(16 * 50),
    w = 1:16,
    y = unlist(hmm_biofam$observations[[1]][1:50, ])
  )
  expect_error(
    fit <- estimate_mnhmm(
      "y", n_states = 4, n_clusters = 2, 
      initial_formula = ~ z, cluster_formula = ~ z, 
      transition_formula = ~w, emission_formula = ~ w, 
      data = d, time = "time", id = "group", maxeval = 1,
      method = "DNM"
    ),
    NA
  )
  expect_error(
    fit <- bootstrap_coefs(fit, nsim = 10),
    NA
  )
  expect_error(
    p <- get_cluster_probs(fit),
    NA
  )
  expect_error(
    p <- get_initial_probs(fit),
    NA
  )
  expect_error(
    p <- get_emission_probs(fit),
    NA
  )
  expect_error(
    p <- get_transition_probs(fit),
    NA
  )
  expect_error(
    coef(fit, probs = 0.5),
    NA
  )
})
test_that("'get_probs' works for 'hmm' and 'mhmm", {
  expect_error(
    get_initial_probs(hmm_biofam),
    NA
  )
  expect_error(
    get_transition_probs(hmm_biofam),
    NA
  )
  expect_error(
    get_emission_probs(hmm_biofam),
    NA
  )
  expect_error(
    get_initial_probs(mhmm_biofam),
    NA
  )
  expect_error(
    get_transition_probs(mhmm_biofam),
    NA
  )
  expect_error(
    get_emission_probs(mhmm_biofam),
    NA
  )
  expect_error(
    get_cluster_probs(mhmm_biofam),
    NA
  )
})
test_that("'get_probs' and 'coef' works for 'fanhmm'", {
  data("fanhmm_leaves")
  expect_error(
    cf <- coef(fanhmm_leaves, probs = c(0.1, 0.9)),
    NA
  )
  expect_equal(names(cf), c("initial", "transition", "emission"))
  expect_equal(
    names(cf$transition), 
    c("state_from", "state_to", "coefficient", "estimate", "q10", "q90")
  )
  expect_equal(
    cf$emission$estimate[1:3], 
    c(-1.30118969295931, 0.829786497113022, 0.471403195846288)
  )
  expect_error(
    p <- get_emission_probs(fanhmm_leaves),
    NA
  )
  expect_equal(
    names(p), 
    c("workplace", "father", "state", "leave", "probability")
  )
  expect_error(
    marginals <- get_marginals(fanhmm_leaves, probs = c(0.1, 0.5, 0.9)),
    NA
  )
  expect_equal(
    names(marginals), 
    c("states", "responses", "transitions", "emissions")
  )
  expect_equal(
    names(marginals$states), 
    c("state", "probability", "q10", "q50", "q90")
  )
  expect_true(
    all(marginals$emissions$q10 < marginals$emissions$probability)
  )
})
