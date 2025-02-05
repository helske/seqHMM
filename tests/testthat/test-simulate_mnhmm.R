test_that("simulate_mnhmm, coef and get_probs works", {
  set.seed(1)
  p <- 20
  n <- 10
  d <- data.frame(
    person = rep(1:p, each = n),
    month = rep(1:n, p),
    x = rnorm(n * p),
    z = rnorm(n * p),
    w = rnorm(n * p)
  )
  expect_error(
    sim <- simulate_mnhmm(
      n_sequences = p, 
      sequence_lengths = as.integer(n), 
      n_symbols = c(4L, 3L), 
      n_states = 2, 
      n_clusters = 3, 
      initial_formula = ~1, transition_formula = ~ x, 
      emission_formula = ~ x + z, cluster_formula = ~w, 
      data = d, time = "month", id = "person", 
      coefs = "random", init_sd = 2),
    NA
  )
  expect_equal(
    table(unlist(sim$states)),
    structure(c(
      `Cluster 1: State 1` = 3L, `Cluster 1: State 2` = 7L, 
      `Cluster 2: State 1` = 6L, `Cluster 2: State 2` = 4L, 
      `Cluster 3: State 1` = 83L, `Cluster 3: State 2` = 97L, 
      `*` = 0L, `%` = 0L), dim = 8L, dimnames = structure(list(
        c("Cluster 1: State 1", "Cluster 1: State 2", "Cluster 2: State 1", 
          "Cluster 2: State 2", "Cluster 3: State 1", "Cluster 3: State 2", 
          "*", "%")), names = ""), class = "table")
    
  )
  expect_error(
    fit <- estimate_mnhmm(
      sim$model$obs, n_states = 2, 
      n_clusters = 3, 
      initial_formula = ~1, transition_formula = ~ x, 
      emission_formula = ~ x + z, cluster_formula = ~w, 
      data = d, time = "month", id = "person", 
      inits = sim$model$etas, maxeval = 1,
      maxeval_em_dnm = 1),
    NA
  )
  expect_error(
    cf <- coef(fit),
    NA
  )
  expect_error(
    p <- get_probs(fit),
    NA
  )
})

