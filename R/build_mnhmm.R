#' Build a Mixture Non-homogeneous Hidden Markov Model
#' @noRd
build_mnhmm <- function(
    observations, n_states, n_clusters, initial_formula, 
    transition_formula, emission_formula, cluster_formula,
    data, time, id, state_names = NULL, channel_names = NULL, 
    cluster_names = NULL) {
  
  stopifnot_(
    !missing(n_clusters) && checkmate::test_int(x = n_clusters, lower = 2L), 
    "Argument {.arg n_clusters} must be a single positive integer larger than 1."
  )
  n_clusters <- as.integer(n_clusters)
  if (is.null(cluster_names)) {
    cluster_names <- paste("Cluster", seq_len(n_clusters))
  } else {
    stopifnot_(
      length(cluster_names) == n_clusters,
      "Length of {.arg cluster_names} is not equal to the number of clusters."
    )
  }
  out <- create_base_nhmm(
    observations, data, time, id, n_states, state_names, channel_names, 
    initial_formula, transition_formula, emission_formula, cluster_formula, 
    cluster_names)
  out$model$call <- match.call()
  structure(
    out$model,
    class = "mnhmm",
    nobs = attr(out$observations, "nobs"),
    df = out$extras$n_pars,
    type = paste0(out$extras$multichannel, "mnhmm_", out$extras$model_type)
  )
}
