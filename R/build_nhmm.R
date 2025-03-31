#' Build a Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_nhmm <- function(
    responses, n_states, initial_formula, transition_formula, 
    emission_formula, data, id_var, time_var, state_names = NULL,
    scale = TRUE) {
  
  out <- create_base_nhmm(
    responses, data, id_var, time_var, n_states, state_names,
    initial_formula, transition_formula, emission_formula, scale = scale) 
  out[c("cluster_names", "n_clusters", "X_omega")] <- NULL
  out$model$etas <- stats::setNames(
    create_initial_values(list(), out$model, 0), c("pi", "A", "B")
  )
  structure(
    out$model,
    class = "nhmm",
    nobs = out$extras$n_obs,
    df = out$extras$np_pi + out$extras$np_A + out$extras$np_B,
    type = paste0(out$extras$multichannel, "nhmm"),
    intercept_only = out$extras$intercept_only,
    np_pi = out$extras$np_pi,
    np_A = out$extras$np_A,
    np_B = out$extras$np_B
  )
}
