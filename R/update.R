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
  X <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable
  )
  object$X_initial <- X$X
  attr(object, "iv_pi") <- X$iv
  X <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths
  )
  object$X_transition <- X$X
  attr(object, "iv_A") <- X$iv
  attr(object, "tv_A") <- X$tv
  X <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$n_channels, object$time_variable, object$id_variable, 
    object$sequence_lengths
  )
  object$X_emission <- X$X
  attr(object, "iv_B") <- X$iv
  attr(object, "tv_B") <- X$tv
  object
}
#' @rdname update_nhmm
#' @export
update.mnhmm <- function(object, newdata, ...) {
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  object$n_sequences <- length(unique(newdata[[object$id_variable]]))
  object$length_of_sequences <- length(unique(newdata[[object$time_variable]]))
  if (!is.null(object$data)) object$data <- newdata
  X <- model_matrix_initial_formula(
    object$initial_formula, newdata, object$n_sequences,
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable
  )
  object$X_initial <- X$X
  attr(object, "iv_pi") <- X$iv
  X <- model_matrix_transition_formula(
    object$transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$time_variable, 
    object$id_variable, object$sequence_lengths
  )
  object$X_transition <- X$X
  attr(object, "iv_A") <- X$iv
  attr(object, "tv_A") <- X$tv
  X <- model_matrix_emission_formula(
    object$emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    object$n_channels, object$time_variable, object$id_variable, 
    object$sequence_lengths
  )
  object$X_emission <- X$X
  attr(object, "iv_B") <- X$iv
  attr(object, "tv_B") <- X$tv
  X <- model_matrix_cluster_formula(
    object$cluster_formula, newdata, object$n_sequences, object$n_clusters,
    object$time_variable, object$id_variable
  )
  object$X_cluster <- X$X
  attr(object, "iv_omega") <- X$iv
  object
}