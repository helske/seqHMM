#' Build a Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_nhmm <- function(
    n_states, emission_formula, initial_formula, transition_formula, 
    data, id_var, time_var, state_names = NULL,
    scale = TRUE) {
  
  out <- create_base_nhmm(
    data, id_var, time_var, n_states, state_names,
    emission_formula, initial_formula, transition_formula, scale = scale) 
  out[c("cluster_names", "n_clusters", "X_omega")] <- NULL
  out$model$etas <- stats::setNames(
    create_initial_values(list(), out$model, 0), c("pi", "A", "B")
  )
  out$model$gammas$pi <- eta_to_gamma_mat(out$model$etas$pi)
  out$model$gammas$A <- eta_to_gamma_cube(out$model$etas$A)
  out$model$gammas$B <- eta_to_gamma_cube_field(out$model$etas$B)
  structure(
    out$model,
    class = "nhmm",
    nobs = out$extras$n_obs,
    df = out$extras$np_pi + out$extras$np_A + out$extras$np_B,
    np_pi = out$extras$np_pi,
    np_A = out$extras$np_A,
    np_B = out$extras$np_B
  )
}
