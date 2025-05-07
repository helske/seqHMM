test_that("simulate_mnhmm and coef works", {
  set.seed(1)
  p <- 50
  n <- 10
  d <- data.frame(
    person = rep(1:p, each = n),
    month = rep(1:n, p),
    x = stats::rnorm(n * p),
    z = stats::rnorm(n * p),
    w = stats::rnorm(n * p),
    y = factor(NA, levels = letters[1:4]),
    g = factor(NA, levels = 1:3)
  )
  expect_error(
    sim <- simulate_mnhmm(
      n_states = 2, 
      n_clusters = 3, 
      initial_formula = ~1, transition_formula = ~ x, 
      emission_formula = c(y, g) ~ x + z, cluster_formula = ~w, 
      data = d, time = "month", id = "person"),
    NA
  )
  expect_equal(
    c(table(sim$states$state)),
    c(`Cluster 1: State 1` = 88L, `Cluster 1: State 2` = 142L, `Cluster 2: State 1` = 140L, 
      `Cluster 2: State 2` = 40L, `Cluster 3: State 1` = 46L, `Cluster 3: State 2` = 44L
    )
  )
  expect_error(
    fit <- estimate_mnhmm(
      n_states = 2, 
      n_clusters = 3, 
      initial_formula = ~1, transition_formula = ~ x, 
      emission_formula = c(y, g) ~ x + z, cluster_formula = ~w, 
      data = sim$data, time = "month", id = "person", 
      inits = sim$model$etas, maxeval = 1,
      maxeval_em_dnm = 1),
    NA
  )
  expect_error(
    cf <- coef(fit),
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
})
