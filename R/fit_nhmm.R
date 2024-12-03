#' Estimate a Non-homogeneous Hidden Markov Model
#'
#' @noRd
fit_nhmm <- function(model, inits, init_sd, restarts, lambda, method,
                     bound, save_all_solutions = FALSE, control_restart = list(), 
                     control_mstep = list(), ...) {
  
  stopifnot_(
    checkmate::test_int(x = restarts, lower = 0L), 
    "Argument {.arg restarts} must be a single integer."
  )
  method <- match.arg(method, c("EM-DNM", "DNM", "EM"))
  stopifnot_(
    checkmate::check_number(lambda, lower = 0), 
    "Argument {.arg lambda} must be a single non-negative {.cls numeric} value."
  )
  stopifnot_(
    checkmate::check_number(bound, lower = 0), 
    "Argument {.arg bound} must be a single non-negative {.cls numeric} value."
  )
  control <- utils::modifyList(
    list(
      ftol_abs = 0,
      xtol_abs = 0,
      ftol_rel = 1e-10,
      xtol_rel = 1e-6,
      maxeval = 1e4,
      print_level = 0,
      algorithm = "NLOPT_LD_LBFGS",
      maxeval_em_dnm = 10
    ),
    list(...)
  )
  control_restart <- utils::modifyList(control, control_restart)
  stopifnot_(
    identical(control$algorithm, control_restart$algorithm),
    c("Cannot mix different algorithms for multistart and final optimization.",
      "Found algorithm {.val {control$algorithm}} for final optimization and 
      {.val {control_restart$algorithm}} for multistart.")
  )
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
    model$gammas$pi <- eta_to_gamma_mat(model$etas$pi)
    model$gammas$A <- eta_to_gamma_cube(model$etas$A)
    if (model$n_channels == 1L) {
      model$gammas$B <- eta_to_gamma_cube(model$etas$B)
    } else {
      model$gammas$B <- eta_to_gamma_cube_field(model$etas$B)
    }
    return(model)
  }
  if (method == "EM-DNM") {
    out <- em_dnm_nhmm(
      model, inits, init_sd, restarts, lambda, bound, control, 
      control_restart, control_mstep, save_all_solutions
    )
  }
  if (method == "DNM") {
    out <- dnm_nhmm(
      model, inits, init_sd, restarts, lambda,  bound, control, control_restart, 
      save_all_solutions 
    )
  }
  if (method == "EM") {
    out <- em_nhmm(
      model, inits, init_sd, restarts, lambda,  bound, control, 
      control_restart, control_mstep, save_all_solutions 
    )
  }
  out
}
