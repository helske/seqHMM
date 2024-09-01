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
  out[c("cluster_names", "n_clusters", "X_cluster", "coef_names_cluster")] <- NULL
  out$model$call <- match.call()
  structure(
    out$model,
    class = "nhmm",
    time_variable = time,
    id_variable = id,
    nobs = attr(out$observations, "nobs"),
    df = out$extras$n_pars,
    type = paste0(out$extras$multichannel, "nhmm_", out$extras$model_type)
  )
}
