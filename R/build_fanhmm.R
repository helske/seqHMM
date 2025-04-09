#' #' Build a Feedback-Augmented Non-homogeneous Hidden Markov Model
#' #'
#' #' @noRd
#' build_fanhmm <- function(
#'     n_states, initial_formula, 
#'     transition_formula, emission_formula, autoregression_formula, 
#'     feedback_formula, data, id_var, time_var, state_names = NULL, 
#'     scale = TRUE) {
#'  
#'   stopifnot_(
#'     is.null(autoregression_formula) || inherits(autoregression_formula, "formula"), 
#'     "Argument {.arg autoregression_formula} must be {.val NULL} or a {.cls formula} object."
#'   )
#'   stopifnot_(
#'     is.null(feedback_formula) || inherits(feedback_formula, "formula"), 
#'     "Argument {.arg feedback_formula} must be {.val NULL} or a {.cls formula} object."
#'   )
#'   stopifnot_(
#'     !is.null(autoregression_formula) || !is.null(feedback_formula),
#'     "Provide {.arg autoregression_formula} and/or {.arg feedback_formula} for FAN-HMM."
#'   )
#'   
#'   if (inherits(emission_formula, "formula")) {
#'     responses <- get_responses(emission_formula)
#'     C <- length(responses)
#'     if (C > 1L) {
#'       rhs <- deparse1(emission_formula[[3L]])
#'       emission_formula <- lapply(
#'         responses, \(y) as.formula(paste(y, " ~ ", rhs), env = environment(emission_formula))
#'       )
#'     } else {
#'       emission_formula <- list(emission_formula)
#'     }
#'   } else {
#'     responses <- vapply(emission_formula, get_responses, allow_mv = FALSE, "")
#'     C <- length(responses)
#'   }
#'   
#'   lag_obs <- paste0("lag_", responses)
#'   if (!is.null(feedback_formula)) {
#'     terms_feedback <- attr(stats::terms(feedback_formula), "term.labels")
#'     if (length(terms_feedback) == 0) {
#'       terms_feedback <- lag_obs
#'     } else {
#'       terms_feedback <- paste(
#'         c(lag_obs, 
#'           paste(
#'             lag_obs, 
#'             terms_feedback, 
#'             sep = " : "
#'           )),
#'         collapse = " + "
#'       )
#'     }
#'     transition_formula <- update(
#'       transition_formula, 
#'       paste("~ . + ", terms_feedback)
#'     )
#'   }
#'   if (ar <- !is.null(autoregression_formula)) {
#'     terms_autoregression <- attr(
#'       stats::terms(autoregression_formula), "term.labels"
#'     )
#'     if (length(terms_autoregression) == 0) {
#'       terms_autoregression <- lag_obs
#'     } else {
#'       terms_autoregression <- paste(
#'         c(lag_obs,
#'           paste(
#'             lag_obs, 
#'             terms_autoregression, 
#'             sep = " : "
#'           )
#'         ), collapse = " + "
#'       )
#'     }
#'     emission_formula <- update(
#'       emission_formula, 
#'       paste("~ . + ", terms_autoregression)
#'     )
#'   }
#' 
#'   out <- create_base_nhmm(
#'     data, id_var, time_var, n_states, state_names, 
#'     emission_formula, initial_formula, transition_formula, scale = scale, 
#'     fanhmm = TRUE)
#'   
#'   out[c("cluster_names", "n_clusters", "X_omega")] <- NULL
#'   out$model$etas <- stats::setNames(
#'     create_initial_values(list(), out$model, 0), c("pi", "A", "B")
#'   )
#'   out$model$gammas$pi <- eta_to_gamma_mat(out$model$etas$pi)
#'   out$model$gammas$A <- eta_to_gamma_cube(out$model$etas$A)
#'   out$model$gammas$B <- eta_to_gamma_cube_field(out$model$etas$B)
#'   structure(
#'     c(
#'       out$model,
#'       list(
#'         autoregression_formula = autoregression_formula, 
#'         feedback_formula = feedback_formula
#'       )
#'     ),
#'     class = c("fanhmm", "nhmm"),
#'     nobs = out$extras$n_obs,
#'     df = out$extras$np_pi + out$extras$np_A + out$extras$np_B,
#'     np_pi = out$extras$np_pi,
#'     np_A = out$extras$np_A,
#'     np_B = out$extras$np_B
#'   )
#' }
