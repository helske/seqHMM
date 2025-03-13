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
  
  if (!is.null(object$data)) object$data <- newdata
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  ids <- unique(newdata[[object$id_variable]])
  object$n_sequences <- length(ids)
  times <- unique(newdata[[object$time_variable]])
  object$length_of_sequences <- length(times)
  object$sequence_lengths <- as.integer(table(newdata[[object$id_variable]]))
  newdata <- fill_time(newdata, object$time_variable, object$id_variable)
  
  if (inherits(object, "fanhmm") & !is.null(object$autoregression_formula)) {
    newdata[[paste0("lag_", object$channel_names)]] <- 
      group_lag(newdata, object$id_variable, object$channel_names)
  }
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
    scale = TRUE, R_inv = attr(object$X_B, "R_inv"), 
    fanhmm = inherits(object, "fanhmm")
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
  #TODO, how to update sequence_lengths allowing NA values in predict?
  # object$sequence_lengths <- pmax(
  #   object$sequence_lengths,
  #   attr(object$observations, "sequence_lengths")
  # )
  object
}
#' @rdname update_nhmm
#' @export
update.mnhmm <- function(object, newdata, ...) {
  if (!is.null(object$data)) object$data <- newdata
  newdata <- .check_data(newdata, object$time_variable, object$id_variable)
  ids <- unique(newdata[[object$id_variable]])
  object$n_sequences <- length(ids)
  times <- unique(newdata[[object$time_variable]])
  object$length_of_sequences <- length(times)
  object$sequence_lengths <- as.integer(table(newdata[[object$id_variable]]))
  newdata <- fill_time(newdata, object$time_variable, object$id_variable)
  
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

update_W_for_fanhmm <- function(object, newdata, ...) {
  
  newdata <- newdata[order(newdata[[object$id_variable]], newdata[[object$time_variable]]), ]
  newdata <- fill_time(newdata, object$time_variable, object$id_variable)
  W_A <- W_B <- vector("list", object$n_symbols)
  newdata[[object$channel_names]] <- factor(object$symbol_names[1], levels = object$symbol_names)
  newdata[[paste0("lag_", object$channel_names)]] <- newdata[[object$channel_names]]
  for (i in seq_len(object$n_symbols)) {
    newdata[[object$channel_names]][] <- object$symbol_names[i]
    newdata[[paste0("lag_", object$channel_names)]][] <- object$symbol_names[i]
    
    W_A[[i]] <- model_matrix_transition_formula(
      object$transition_formula, newdata, object$n_sequences, 
      object$length_of_sequences, object$n_states, object$time_variable, 
      object$id_variable, object$sequence_lengths,
      X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
      scale = TRUE, R_inv = attr(object$X_A, "R_inv")
    )$X
    W_B[[i]] <- model_matrix_emission_formula(
      object$emission_formula, newdata, object$n_sequences, 
      object$length_of_sequences, object$n_states, object$n_symbols, 
      object$time_variable, object$id_variable, 
      object$sequence_lengths,
      X_mean = attr(object$X_B, "X_mean"), check = FALSE, 
      scale = TRUE, R_inv = attr(object$X_B, "R_inv"), 
      fanhmm = inherits(object, "fanhmm")
    )$X
  }
  list(W_A = W_A, W_B = W_B)
}
