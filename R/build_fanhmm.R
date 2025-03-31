#' Build a Feedback-Augmented Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_fanhmm <- function(
    responses, n_states, initial_formula, 
    transition_formula, emission_formula, autoregression_formula, 
    feedback_formula, data, id_var, time_var, state_names = NULL, scale = TRUE) {
  
 
  stopifnot_(
    !missing(responses) && checkmate::test_character(x = responses), 
    "Argument {.arg responses} must be a character vector defining the response 
    variable(s) in the {.arg data}."
  )
  stopifnot_(
    length(responses) == 1L, 
    "Currently only single-channel responses are supported for FAN-HMM."
  )
  stopifnot_(
    is.null(autoregression_formula) || inherits(autoregression_formula, "formula"), 
    "Argument {.arg autoregression_formula} must be {.val NULL} or a {.cls formula} object."
  )
  stopifnot_(
    is.null(feedback_formula) || inherits(feedback_formula, "formula"), 
    "Argument {.arg feedback_formula} must be {.val NULL} or a {.cls formula} object."
  )
  stopifnot_(
    !is.null(autoregression_formula) || !is.null(feedback_formula),
    "Provide {.arg autoregression_formula} and/or {.arg feedback_formula} for FAN-HMM."
  )
  
  lag_obs <- paste0("lag_", responses)
  if (!is.null(autoregression_formula)) {
    terms_autoregression <- attr(
      stats::terms(autoregression_formula), "term.labels"
    )
    if (length(terms_autoregression) == 0) {
      terms_autoregression <- lag_obs
    } else {
      terms_autoregression <- paste(
        c(lag_obs,
          paste(
            lag_obs, 
            terms_autoregression, 
            sep = ":"
          )
        ), collapse = "+"
      )
    }
    emission_formula <- update(
      emission_formula, 
      paste("~ . + ", terms_autoregression)
    )
  }
  if (!is.null(feedback_formula)) {
    terms_feedback <- attr(stats::terms(feedback_formula), "term.labels")
    if (length(terms_feedback) == 0) {
      terms_feedback <- lag_obs
    } else {
      terms_feedback <- paste(
        c(lag_obs, 
          paste(
            lag_obs, 
            terms_feedback, 
            sep = ":"
          )),
        collapse = "+"
      )
    }
    transition_formula <- update(
      transition_formula, 
      paste("~ . + ", terms_feedback)
    )
  }
  out <- create_base_nhmm(
    responses, data, id_var, time_var, n_states, state_names, 
    initial_formula, transition_formula, emission_formula, scale = scale, 
    fanhmm = TRUE)
  
  out[c("cluster_names", "n_clusters", "X_omega")] <- NULL
  out$model$etas <- stats::setNames(
    create_initial_values(list(), out$model, 0), c("pi", "A", "B")
  )
  structure(
    c(
      out$model,
      list(
        autoregression_formula = autoregression_formula, 
        feedback_formula = feedback_formula
      )
    ),
    class = c("fanhmm", "nhmm"),
    nobs = out$extras$n_obs,
    df = out$extras$np_pi + out$extras$np_A + out$extras$np_B,
    type = paste0(out$extras$multichannel, "fanhmm"),
    intercept_only = out$extras$intercept_only,
    np_pi = out$extras$np_pi,
    np_A = out$extras$np_A,
    np_B = out$extras$np_B
  )
}
