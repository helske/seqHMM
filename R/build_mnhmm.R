#' Build a Mixture Non-homogeneous Hidden Markov Model
#' @noRd
build_mnhmm <- function(
    responses, n_states, n_clusters, initial_formula, 
    transition_formula, emission_formula, cluster_formula,
    data, id_var, time_var, state_names = NULL, cluster_names = NULL, scale = TRUE) {
  
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
    responses, data, id_var, time_var, n_states, state_names, 
    initial_formula, transition_formula, emission_formula, cluster_formula, 
    cluster_names, scale = scale)
  out$model$etas <- stats::setNames(
    create_initial_values(list(), out$model, 0), c("pi", "A", "B", "omega")
  )
  structure(
    out$model,
    class = "mnhmm",
    nobs = out$extras$n_obs,
    df = out$extras$np_omega + 
      n_clusters * (out$extras$np_pi + out$extras$np_A + out$extras$np_B),
    type = paste0(out$extras$multichannel, "mnhmm"),
    intercept_only = out$extras$intercept_only,
    np_pi = n_clusters * out$extras$np_pi,
    np_A = n_clusters * out$extras$np_A,
    np_B = n_clusters * out$extras$np_B,
    np_omega = out$extras$np_omega
  )
}
