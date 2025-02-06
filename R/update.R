#' Update Covariate Values of NHMM
#' 
#' This function can be used to replace original covariate values of NHMMs. 
#' The responses, model formulae and estimated coefficients are not altered.
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param newdata A data frame containing the new covariate values.
#' @param ... Ignored.
#' @rdname update_nhmm
#' @export
update.nhmm <- function(object, newdata, ...) {
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  if (!is.null(object$data)) object$data <- newdata
  
  ids <- unique(newdata[[object$id_variable]])
  object$n_sequences <- length(ids)
  times <- unique(newdata[[object$time_variable]])
  object$length_of_sequences <- length(times)
  if (inherits(object, "fanhmm") & !is.null(object$autoregression_formula)) {
    newdata[[paste0("lag_", object$channel_names)]] <- 
      group_lag(newdata, object$id_variable, object$channel_names)
  }
  object$X_pi <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, 
    X_mean = attr(object$X_pi, "X_mean"), check = FALSE, 
    scale = TRUE, R = attr(object$X_pi, "R")
  )$X
  object$X_A <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths,
    X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
    scale = TRUE, R = attr(object$X_A, "R")
  )$X
  object$X_B <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$time_variable, object$id_variable, 
    object$sequence_lengths,
    X_mean = attr(object$X_B, "X_mean"), check = FALSE, 
    scale = TRUE, R = attr(object$X_B, "R"), fanhmm = inherits(object, "fanhmm")
  )$X
  observations <- lapply(
    object$channel_names,
    function(y) {
      stopifnot_(
        !is.null(newdata[[y]]), 
        "Can't find response variable {.var {y}} in {.arg newdata}."
      )
      stopifnot_(
        is.factor(newdata[[y]]), 
        "Response {.var {y}} in {.arg newdata} should be a factor."
      )
      x <- suppressMessages(
        seqdef(matrix(
          newdata[[y]], 
          object$n_sequences, 
          object$length_of_sequences, byrow = TRUE),
          id = ids,
          alphabet = alphabet(object$observations)
        )
      )
      colnames(x) <- sort(times)
      x
    }
  )
  object$observations <- .check_observations(observations, object$channel_names)
  object$sequence_lengths <- attr(object$observations, "sequence_lengths")
  object
}
#' @rdname update_nhmm
#' @export
update.mnhmm <- function(object, newdata, ...) {
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  if (!is.null(object$data)) object$data <- newdata
  
  ids <- unique(newdata[[object$id_variable]])
  object$n_sequences <- length(ids)
  times <- unique(newdata[[object$time_variable]])
  object$length_of_sequences <- length(times)
  
  object$X_pi <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable,
    X_mean = attr(object$X_pi, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_pi, "R_inv")
  )$X
  object$X_A <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths,
    X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_A, "R_inv")
  )$X
  object$X_B <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$time_variable, object$id_variable, 
    object$sequence_lengths,
    X_mean = attr(object$X_B, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_B, "R_inv")
  )$X
  object$X_omega <- model_matrix_cluster_formula(
    object$cluster_formula, newdata, object$n_sequences, object$n_clusters,
    object$time_variable, object$id_variable,
    X_mean = attr(object$X_omega, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_omega, "R_inv")
  )$X
  observations <- lapply(
    object$channel_names,
    function(y) {
      stopifnot_(
        !is.null(newdata[[y]]), 
        "Can't find response variable {.var {y}} in {.arg newdata}."
      )
      stopifnot_(
        is.factor(newdata[[y]]), 
        "Response {.var {y}} in {.arg newdata} should be a factor."
      )
      x <- suppressMessages(
        seqdef(matrix(
          newdata[[y]], 
          object$n_sequences, 
          object$length_of_sequences, byrow = TRUE),
          id = ids,
          alphabet = alphabet(object$observations)
        )
      )
      colnames(x) <- sort(times)
      x
    }
  )
  object$observations <- .check_observations(observations, object$channel_names)
  object
}
