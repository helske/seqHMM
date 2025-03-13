#' Build a Feedback-Augmented Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_fanhmm <- function(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, autoregression_formula, 
    feedback_formula, data, time, id, state_names = NULL, scale = TRUE) {
  
  y_in_data <- checkmate::test_character(observations)
  stopifnot_(
    y_in_data && !is.null(data[[observations]]),
    "For FAN-HMM, the response variable {.arg observations} must be in the {.arg data}."
  )
  stopifnot_(
    length(observations) == 1L,
    "Currently only single-channel responses are supported for FAN-HMM.")
  
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
  stopifnot_(
    inherits(initial_formula, "formula"), 
    "Argument {.arg initial_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(transition_formula, "formula"), 
    "Argument {.arg transition_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(emission_formula, "formula"), 
    "Argument {.arg emission_formula} must be a {.cls formula} object.")
  
  data <- .check_data(data, time, id)
  data <- fill_time(data, time, id)
  if (!is.null(autoregression_formula)) {
    terms_autoregression <- attr(
      stats::terms(autoregression_formula), "term.labels"
    )
    if (length(terms_autoregression) == 0) {
      terms_autoregression <- paste0("lag_", observations)
    } else {
      terms_autoregression <- paste(
        c(paste0("lag_", observations),
          paste(
            paste0("lag_", observations), 
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
    data[[paste0("lag_", observations)]] <- group_lag(data, id, observations)
  }
  if (!is.null(feedback_formula)) {
    terms_feedback <- attr(stats::terms(feedback_formula), "term.labels")
    if (length(terms_feedback) == 0) {
      terms_feedback <- observations
    } else {
      terms_feedback <- paste(
        c(observations, 
          paste(
            observations, 
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
    observations, data, time, id, n_states, state_names, channel_names = NULL,
    initial_formula, transition_formula, emission_formula, scale = scale, 
    check_formulas = FALSE, fanhmm = !is.null(autoregression_formula))
  stopifnot_(
    !any(out$model$observations == attr(out$model$observations, "nr")),
    "FAN-HMM does not yet support missing values in the observations."
  )
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
    nobs = attr(out$model$observations, "nobs") - out$model$n_sequences,
    df = out$extras$np_pi + out$extras$np_A + out$extras$np_B,
    type = paste0(out$extras$multichannel, "fanhmm"),
    intercept_only = out$extras$intercept_only,
    np_pi = out$extras$np_pi,
    np_A = out$extras$np_A,
    np_B = out$extras$np_B
  )
}
