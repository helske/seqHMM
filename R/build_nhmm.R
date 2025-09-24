#' Build a Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_nhmm <- function(
    n_states, emission_formula, initial_formula, transition_formula, 
    data, id_var, time_var, state_names = NULL, scale = TRUE,
    prior_obs = "fixed", coefs = list(), check = NULL, drop_levels = TRUE) {
  
  model <- create_base_nhmm(
    data, id_var, time_var, n_states, state_names,
    emission_formula, initial_formula, transition_formula, scale = scale,
    prior_obs = prior_obs, check = check, drop_levels = drop_levels) 
  model[c("cluster_names", "X_omega")] <- NULL
  model$etas <- create_initial_values(coefs, model, 0)
  model$gammas$gamma_pi <- eta_to_gamma_mat(model$etas$eta_pi)
  model$gammas$gamma_A <- eta_to_gamma_cube(model$etas$eta_A)
  model$gammas$gamma_B <- drop(eta_to_gamma_cube_field(model$etas$eta_B))
  model
}
