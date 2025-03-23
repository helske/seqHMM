test_that("simulate_mnhmm and coef works", {
  set.seed(1)
  p <- 20
  n <- 10
  d <- data.frame(
    person = rep(1:p, each = n),
    month = rep(1:n, p),
    x = stats::rnorm(n * p),
    z = stats::rnorm(n * p),
    w = stats::rnorm(n * p)
  )
  expect_error(
    sim <- simulate_mnhmm(
      responses = c("y", "g"),
      n_sequences = p, 
      sequence_lengths = as.integer(n), 
      n_symbols = c(4L, 3L), 
      n_states = 2, 
      n_clusters = 3, 
      initial_formula = ~1, transition_formula = ~ x, 
      emission_formula = ~ x + z, cluster_formula = ~w, 
      data = d, time = "month", id = "person"),
    NA
  )
  expect_equal(
    c(table(sim$states$state)),
    c(`Cluster 1: State 1` = 2L, `Cluster 1: State 2` = 8L, `Cluster 2: State 1` = 5L, 
      `Cluster 2: State 2` = 5L, `Cluster 3: State 1` = 92L, `Cluster 3: State 2` = 88L
    )
  )
  expect_error(
    fit <- estimate_mnhmm(
      c("y", "g"), n_states = 2, 
      n_clusters = 3, 
      initial_formula = ~1, transition_formula = ~ x, 
      emission_formula = ~ x + z, cluster_formula = ~w, 
      data = sim$model$data, time = "month", id = "person", 
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
