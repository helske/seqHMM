set.seed(1)
p <- 100
n <- 20
d <- data.frame(
  person = rep(1:p, each = n),
  month = rep(1:n, p),
  x = rnorm(n * p),
  z = rnorm(n * p),
  w = rnorm(n * p)
)
test_that("simulate_mnhmm works", {
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
})

