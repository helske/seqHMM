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
  object$X_pi <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, 
    attr(object$X_pi, "X_mean"), attr(object$X_pi, "X_sd"), FALSE
  )$X
  object$X_A <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths,
    attr(object$X_A, "X_mean"), attr(object$X_A, "X_sd"), FALSE
  )$X
  object$X_B <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$n_channels, object$time_variable, object$id_variable, 
    object$sequence_lengths,
    attr(object$X_B, "X_mean"), attr(object$X_B, "X_sd"), FALSE
  )$X
  object
}
#' @rdname update_nhmm
#' @export
update.mnhmm <- function(object, newdata, ...) {
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  if (!is.null(object$data)) object$data <- newdata
  object$X_pi <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable,
    attr(object$X_pi, "X_mean"), attr(object$X_pi, "X_sd"), FALSE
  )$X
  object$X_A <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths,
    attr(object$X_A, "X_mean"), attr(object$X_A, "X_sd"), FALSE
  )$X
  object$X_B <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$n_channels, object$time_variable, object$id_variable, 
    object$sequence_lengths,
    attr(object$X_B, "X_mean"), attr(object$X_B, "X_sd"), FALSE
  )$X
  object$X_omega <- model_matrix_cluster_formula(
    object$cluster_formula, newdata, object$n_sequences, object$n_clusters,
    object$time_variable, object$id_variable,
    attr(object$X_omega, "X_mean"), attr(object$X_omega, "X_sd"), FALSE
  )$X
  object
}
