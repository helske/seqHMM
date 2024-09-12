run_extended_tests <- identical(Sys.getenv("SEQHMM_EXTENDED_TESTS"), "true")

test_that("build_mhmm and estimate_mnhmm give comparable results", {
  skip_if_not(run_extended_tests)
  data("mvad", package = "TraMineR")
  
  d <- reshape(mvad, direction = "long", varying = list(15:86), 
               v.names = "activity", timevar = "month", idvar = "person")

  fit1 <- fit_model(
    build_mhmm(seqdef(mvad[, 15:86]),
               n_states = c(3, 3),
               data = d |> dplyr::filter(time == 1), 
               formula = ~gcse5eq
    ), control_em = list(restart = list(times = 1000))
  )
  fit2 <- estimate_mnhmm(
    "activity", n_states = 3, n_clusters = 2,
    data = d, time = "time", id = "id", 
    cluster_formula = ~ gcse5eq,
    inits = "random", restart = 1000,
  )
  
  expect_equal(logLik(fit1), logLik(fit2), tolerance = 0.1)
  
})
