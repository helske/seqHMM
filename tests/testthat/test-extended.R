run_extended_tests <- identical(Sys.getenv("SEQHMM_EXTENDED_TESTS"), "true")

test_that("build_mhmm and estimate_mnhmm give comparable results", {
  skip_if_not(run_extended_tests)
  data("mvad", package = "TraMineR")
  
  d <- reshape(mvad, direction = "long", varying = list(15:86), 
               v.names = "activity", timevar = "month", idvar = "person")
  
  expect_warning(
    model <- build_mhmm(
      seqdef(mvad[, 15:86]),
      n_states = c(3, 3),
      data = d[d$month == 1, ], 
      formula = ~gcse5eq
    ),
    "Time indices \\(column names\\) of sequences are not coarceable to numeric\\. Replacing them with integers\\."
  )
  fit1 <- fit_model(model, control_em = list(restart = list(times = 10)))
  fit2 <- estimate_mnhmm(
    "activity", n_states = 3, n_clusters = 2, data = d, time = "month", 
    id = "person", cluster_formula = ~ gcse5eq, restarts = 10,
    inits = model
  )
  
  expect_equal(logLik(fit1$model), logLik(fit2), tolerance = 0.1)
  
})

test_that("estimate_nhmm recovers true parameters", {
  skip_if_not(run_extended_tests)
  
  t_seq <- 1e5
  n_seq <- 50
  set.seed(1)
  d <- data.frame(
    person = rep(seq_len(n_seq), each = t_seq),
    year = rep(seq_len(t_seq), n_seq),
    x = stats::rnorm(n_seq * t_seq)
  )
  # Intercepts
  pi <- c(0.7, 0.2, 0.1)
  A <- matrix(c(0.8, 0.1, 0.1,
                0.1, 0.7, 0.2,
                0.1, 0.1, 0.8), S, S, byrow = TRUE)
  B <- matrix(c(0.7, 0.1, 0.1, 0.1,
                0.1, 0.5, 0.2, 0.2,
                0.1, 0.1, 0.4, 0.4), S, M, byrow = TRUE)
  
  # working regression coefficients eta
  coefs <- create_initial_values_(
    inits = list(initial_probs = pi, transition_probs = A, emission_probs = B), 
    init_sd = 0, S, M, K_pi = 1, K_A = 2, K_B = 1
  )
  # transform sum-to-zero constrained gammas to working parameters eta
  tQs <- t(create_Q(S))
  # effect of x in state 1, 2, and 3
  coefs$A[, 2, 1] <- tQs %*% c(1, 0, -1)
  coefs$A[, 2, 2] <- tQs %*% c(0, -1, 1)
  coefs$A[, 2, 3] <- tQs %*% c(-1, 1, 0)
  
  # simulate new sequence data
  sim <- simulate_nhmm(
    n_sequences = n_seq,
    sequence_lengths = t_seq,
    n_symbols = M,
    n_states = S,
    initial_formula = ~ 1,
    transition_formula = ~ x,
    data = d,
    time = "year",
    id = "person",
    coefs = coefs,
    init_sd = 0
  )
  d$y <- c(t(sim$model$observations))
  
  fit <- estimate_nhmm(
    observations = "y",
    n_states = S,
    initial_formula = ~ 1,
    transition_formula = ~ x,
    emission_formula = ~ 1,
    data = d,
    time = "year",
    id = "person",
    lambda = 0,
    method = "DNM",
    inits = coefs, 
    init_sd = 0,
    restarts = 10
  )
  
})
