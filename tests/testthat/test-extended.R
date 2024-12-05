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
