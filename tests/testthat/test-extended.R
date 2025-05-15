run_extended_tests <- identical(Sys.getenv("SEQHMM_EXTENDED_TESTS"), "true")

test_that(
  "parameters are recovered and `mnhmm` and `mhmm` give equivalent results", {
    skip_if_not(run_extended_tests)
    
    set.seed(123)
    n_seq <- 5000
    t_seq <- 50
    pi <- list(c(0.6, 0.4), c(0.8, 0.2))
    omega = c(0.3, 0.7)
    A <- list(diag(0.8, 2) + 0.1, matrix(c(0.95, 0.2, 0.05, 0.8), 2, 2))
    B <- list(
      list(rbind(c(0.5,0.2,0.2,0.1), c(0.1, 0.1, 0.7, 0.1))), 
      list(rbind(c(0.1,0.1,0.4,0.4), c(0.2, 0.3, 0.1, 0.4)))
    )
    expect_error(
      sim <- simulate_mnhmm(
        2, 2, y ~ 1,
        data = data.frame(
          id = rep(1:n_seq, each = t_seq),
          time = rep_len(1:t_seq, n_seq),
          y = factor(sample(letters[1:4], t_seq * n_seq, replace = TRUE), 
                     levels = letters[1:4])),
        id = "id", time = "time",
        coefs = list(cluster_probs = omega, initial_probs = pi, 
                     transition_probs = A, emission_probs = B)
      ),
      NA
    )
    set.seed(1)
    expect_error(
      fit_dnm <- estimate_mnhmm(
        n_states = 2, n_clusters = 2,
        data = sim$data, time = "time", id = "id", 
        emission_formula = y ~ 1, 
        inits = sim$model$etas, init_sd = 0.1,
        method = "DNM", restarts = 50
      ),
      NA
    )
    set.seed(1)
    expect_error(
      fit_em <- estimate_mnhmm(
        n_states = 2, n_clusters = 2,
        data = sim$data, time = "time", id = "id", 
        emission_formula = y ~ 1, 
        inits = sim$model$etas, init_sd = 0.1,
        method = "EM", restarts = 50
      ),
      NA
    )
    
    expect_error(
      init <- build_mhmm(
        observations = data_to_stslist(sim$data, "id", "time", "y"),
        transition_probs = A,
        emission_probs = list(B[[1]][[1]], B[[2]][[1]]),
        initial_probs = pi,
        coefficients = matrix(c(0, sim$model$gammas$gamma_omega[2]*2), 1, 2)
      ),
      NA
    )
    set.seed(1)
    expect_error(
      suppressWarnings(fit_old <- fit_model(
        init, 
        control_em = list(restart = list(times = 50))
      )$model
      )
      ,
      NA
    )
    
    expect_equal(logLik(fit_dnm), logLik(fit_em), tolerance = 0.1)
    expect_equal(logLik(fit_dnm), logLik(fit_old), tolerance = 0.1)
    
    expect_error(coef_true <- coef(sim$model), NA)
    expect_error(coef_est <- coef(fit_dnm), NA)
    expect_equal(coef_est$cluster, coef_true$cluster, tolerance = 0.1)
    expect_equal(coef_est$initial, coef_true$initial, tolerance = 0.1)
    expect_equal(coef_est$emission, coef_true$emission, tolerance = 0.1)
    expect_equal(coef_est$transition, coef_true$transition, tolerance = 0.1)
  }
)

test_that("Mixture FAN-HMM works properly", {
  skip_if_not(run_extended_tests)
  set.seed(123)
  n_seq <- 500
  t_seq <- 20
  pi <- list(c(0.6, 0.4), c(0.8, 0.2))
  omega = c(0.3, 0.7)
  A <- list(diag(0.8, 2) + 0.1, matrix(c(0.95, 0.2, 0.05, 0.8), 2, 2))
  B <- list(
    list(rbind(c(0.5,0.2,0.2,0.1), c(0.1, 0.1, 0.7, 0.1))), 
    list(rbind(c(0.1,0.1,0.4,0.4), c(0.2, 0.3, 0.1, 0.4)))
  )
  expect_error(
    sim <- simulate_mnhmm(
      2, 2, y ~ lag(y),
      data = data.frame(
        id = rep(1:n_seq, each = t_seq),
        time = rep_len(1:t_seq, n_seq),
        y = factor(sample(letters[1:4], t_seq * n_seq, replace = TRUE), 
                   levels = letters[1:4])),
      id = "id", time = "time",
      coefs = list(cluster_probs = omega, initial_probs = pi, 
                   transition_probs = A, emission_probs = B),
      init_sd = 1
    ),
    NA
  )
  y_prior <- prop.table(table(sim$data$y[sim$data$time == 1]))
  set.seed(1)
  expect_error(
    fit_dnm <- estimate_mnhmm(
      n_states = 2, n_clusters = 2,
      data = sim$data, time = "time", id = "id", 
      emission_formula = y ~ lag(y), 
      inits = sim$model$etas, init_sd = 0.1, lambda = 0.1, restarts = 20,
      method = "DNM", prior_obs = list(y = y_prior)
    ),
    NA
  )
  set.seed(1)
  expect_error(
    fit_em <- estimate_mnhmm(
      n_states = 2, n_clusters = 2,
      data = sim$data, time = "time", id = "id", 
      emission_formula = y ~ lag(y), 
      inits = sim$model$etas, init_sd = 0.1, lambda = 0.1, restarts = 20,
      method = "EM", prior_obs = list(y = y_prior)
    ),
    NA
  )
  expect_equal(logLik(fit_dnm), logLik(fit_em),tolerance = 0.1)
  expect_error(coef_dnm <- coef(fit_dnm), NA)
  expect_error(coef_em <- coef(fit_em), NA)
  expect_equal(coef_dnm$cluster, coef_em$cluster, tolerance = 0.1)
  expect_equal(coef_dnm$initial, coef_em$initial, tolerance = 0.1)
  expect_equal(coef_dnm$emission, coef_em$emission, tolerance = 0.1)
  expect_equal(coef_dnm$transition, coef_em$transition, tolerance = 0.1)
})
