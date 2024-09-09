#' Update Covariate Values of NHMM
#' 
#' This function can be used to replace original covariate values of NHMMs. 
#' The model formulae and estimated coefficients are not altered.
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param newdata A data frame containing the new covariate values.
#' @param ... Ignored.
#' @rdname update_nhmm
#' @export
update.nhmm <- function(object, newdata, ...) {
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  object$n_sequences <- length(unique(newdata[[object$id_variable]]))
  object$length_of_sequences <- length(unique(newdata[[object$time_variable]]))
  if (!is.null(object$data)) object$data <- newdata
  object$X_initial <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable
  )$X
  object$X_transition <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths
  )$X
  object$X_emission <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$n_channels, object$time_variable, object$id_variable, 
    object$sequence_lengths
  )$X
  object
}
#' @rdname update_nhmm
#' @export
update.mnhmm <- function(object, newdata, ...) {
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  object$n_sequences <- length(unique(newdata[[object$id_variable]]))
  object$length_of_sequences <- length(unique(newdata[[object$time_variable]]))
  if (!is.null(object$data)) object$data <- newdata
  object$X_initial <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable
  )$X
  object$X_transition <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths
  )$X
  object$X_emission <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$n_channels, object$time_variable, object$id_variable, 
    object$sequence_lengths
  )$X
  object$X_cluster <- model_matrix_cluster_formula(
    object$cluster_formula, newdata, object$n_sequences, object$n_clusters,
    object$time_variable, object$id_variable
  )$X
  object
}