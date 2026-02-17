test_that("SQUAREM works for nhmm", {
  
  set.seed(123)
  n <- 30
  p <- 15
  d <- data.frame(
    person = rep(1:p, each = n),
    time = rep(1:n, p),
    x = rnorm(n * p), 
    z = rnorm(n * p),
    y = factor(NA, levels = letters[4:6])
  )
  
  # Simulate nhmm data with formulas
  sim <- simulate_nhmm(
    n_states = 3,
    initial_formula = ~1,
    transition_formula = ~ x,
    emission_formula = y ~ z,
    data = d,
    id = "person",
    time = "time"
  )
  
  set.seed(386)
  fit_standard <- estimate_nhmm(
    n_states = 3,
    initial_formula = ~1,
    transition_formula = ~ x,
    emission_formula = y ~ z,
    data = sim$data,
    id = "person",
    time = "time",
    method = "EM",
    use_squarem = FALSE
  )
  set.seed(386)
  fit_squarem <- estimate_nhmm(
    n_states = 3,
    initial_formula = ~1,
    transition_formula = ~ x,
    emission_formula = y ~ z,
    data = sim$data,
    id = "person",
    time = "time",
    method = "EM",
    use_squarem = TRUE
  )
  
  expect_true(abs(logLik(fit_standard) - logLik(fit_squarem)) < 1)
  
  expect_true(fit_squarem$estimation_results$iterations <= 
                fit_standard$estimation_results$iterations + 2)
  
  expect_s3_class(fit_standard, "nhmm")
  expect_s3_class(fit_squarem, "nhmm")
})
