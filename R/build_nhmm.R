#' Build a Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_nhmm <- function(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, 
    data, time, id, state_names = NULL, channel_names = NULL) {
  
  out <- create_base_nhmm(
    observations, data, time, id, n_states, state_names, channel_names,
    initial_formula, transition_formula, emission_formula) 
  out[c("cluster_names", "n_clusters", "X_cluster")] <- NULL
  structure(
    out$model,
    class = "nhmm",
    nobs = attr(out$model$observations, "nobs"),
    df = out$extras$np_pi + out$extras$np_A + out$extras$np_B,
    type = paste0(out$extras$multichannel, "nhmm"),
    intercept_only = out$extras$intercept_only,
    np_pi = out$extras$np_pi,
    np_A = out$extras$np_A,
    np_B = out$extras$np_B
  )
}
