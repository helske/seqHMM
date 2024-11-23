#' Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_mnhmm <- function(model, inits, init_sd, restarts, lambda, method, 
                      pseudocount, save_all_solutions = FALSE, 
                      control_restart = list(), control_mstep = list(), ...) {
  
  stopifnot_(
    checkmate::test_int(x = restarts, lower = 0L), 
    "Argument {.arg restarts} must be a single integer."
  )
  control <- utils::modifyList(
    list(
      ftol_abs = 1e-8,
      ftol_rel = 1e-8,
      xtol_abs = 1e-4,
      xtol_rel = 1e-4,
      maxeval = 1e4,
      print_level = 0,
      algorithm = "NLOPT_LD_LBFGS"
    ),
    list(...)
  )
  control_restart <- utils::modifyList(control, control_restart)
  control_mstep <- utils::modifyList(control, control_mstep)
  
  if (identical(inits, "random")) {
    inits <- list(
      initial_probs = NULL, 
      transition_probs = NULL, 
      emission_probs = NULL)
  } else {
    if (is.null(inits$initial_probs)) inits$initial_probs <- NULL
    if (is.null(inits$transition_probs)) inits$transition_probs <- NULL
    if (is.null(inits$emission_probs)) inits$emission_probs <- NULL
  }
  
  if (isTRUE(control$maxeval < 0)) {
    model$etas <- create_initial_values(inits, model, init_sd)
    model$gammas$pi <- c(eta_to_gamma_mat_field(model$etas$pi))
    model$gammas$A <- c(eta_to_gamma_cube_field(model$etas$A))
    if (C == 1L) {
      model$gammas$B <- c(eta_to_gamma_cube_field(model$etas$B))
    } else {
      l <- lengths(model$etas$B)
      gamma_B <- c(eta_to_gamma_cube_field(unlist(model$etas$B, recursive = FALSE)))
      model$gammas$B <- unname(split(gamma_B, rep(seq_along(l), l)))
    }
    model$gammas$omega <- eta_to_gamma_mat(model$etas$omega)
    return(model)
  }
  all_solutions <- NULL
  if (method == "LBFGS") {
    out <- lbfgs_mnhmm(
      model, inits, init_sd, restarts, lambda, control, control_restart, 
      save_all_solutions 
    )
  }
  if (method == "EM") {
    out <- em_mnhmm(
      model, inits, init_sd, restarts, lambda, pseudocount, control, 
      control_restart, control_mstep, save_all_solutions 
    )
  }
  out
}
