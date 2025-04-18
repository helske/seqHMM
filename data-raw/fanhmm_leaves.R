library(seqHMM)
set.seed(1)
library(future)
plan(multisession, workers = 4)

data("leaves")
inits <- list(
  initial_probs = c(0.4, 0.4, 0.2),
  transition_probs = diag(0.85, 3) + 0.05,
  emission_probs = list(leave = matrix(
    c(0.1, 0.2, 0.7,
      0.1, 0.6, 0.3, 
      0.4, 0.4, 0.2
    ), 
    3, 3, byrow = TRUE)
  )
)
progressr::with_progress(
  fanhmm_leaves <- estimate_nhmm(
    n_states = 3,
    emission_formula = 
      leave ~ reform2013 * occupation + lag(leave) * same_occupation,
    initial_formula = ~ reform2013,
    transition_formula =
      ~ lag(leave) * lag_occupation + lag(leave) * lag_reform2013,
    data = leaves, 
    time = "father", id = "workplace", 
    lambda = 0.1, method = "DNM",
    inits = inits,
    init_sd = 2,
    restarts = 50
  )
)

progressr::with_progress(
  fanhmm_leaves <- bootstrap_coefs(fanhmm_leaves, nsim = 100)
)
usethis::use_data(fanhmm_leaves)
