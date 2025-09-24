#' Build a Mixture Non-homogeneous Hidden Markov Model
#' @noRd
build_mnhmm <- function(
    n_states, n_clusters, emission_formula, initial_formula, 
    transition_formula, cluster_formula, data, id_var, time_var, 
    state_names = NULL, cluster_names = NULL, scale = TRUE,
    prior_obs = "fixed", coefs = list(), check = NULL, drop_levels = TRUE) {
  
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
  cluster_names <- factor(cluster_names)
  model <- create_base_nhmm(
    data, id_var, time_var, n_states, state_names, 
    emission_formula, initial_formula, transition_formula, cluster_formula, 
    cluster_names, scale = scale, prior_obs = prior_obs, check = check,
    drop_levels = drop_levels)
  model$etas <- create_initial_values(coefs, model, 0)
  model$gammas$gamma_pi <- drop(eta_to_gamma_mat_field(model$etas$eta_pi))
  model$gammas$gamma_A <- drop(eta_to_gamma_cube_field(model$etas$eta_A))
  model$gammas$gamma_B <- split(
    eta_to_gamma_cube_2d_field(model$etas$eta_B), seq_len(n_clusters)
  )
  model$gammas$gamma_omega <- eta_to_gamma_mat(model$etas$eta_omega)
  model
}
