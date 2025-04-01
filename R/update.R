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
  # avoid CRAN check warning due to NSE
  tmax <- y <- NULL
  initial_formula <- object$initial_formula
  attr(initial_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(initial_formula), object$data
  )
  transition_formula <- object$transition_formula
  attr(transition_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(transition_formula), object$data
  )
  emission_formula <- object$emission_formula
  attr(emission_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(emission_formula), object$data
  )
  time_var <- object$time_variable
  id_var <- object$id_variable
  responses <- object$responses
  newdata <- .check_data(newdata, id_var, time_var, responses)
  ids <- unique(newdata[[id_var]])
  object$n_sequences <- length(ids)
  times <- unique(newdata[[time_var]])
  object$length_of_sequences <- length(times)
  n_obs <- sum(!is.na(newdata[, y, env = list(y = I(responses))])) / object$n_channels
  attr(object, "nobs") <- n_obs
  newdata[, tmax := max(time_var), by = id_var, 
          env = list(id_var = id_var, time_var = time_var)]
  newdata <- fill_time(newdata, id_var, time_var)
  object$sequence_lengths <- newdata[, 
                                     sum(time_var <= tmax, na.rm = TRUE), 
                                     by = id_var, 
                                     env = list(
                                       id_var = id_var, 
                                       time_var = time_var
                                     )]$V1
  newdata[, tmax := NULL]
  if (inherits(object, "fanhmm") && !is.null(object$autoregression_formula)) {
    newdata[[paste0("lag_", responses)]] <- 
      group_lag(newdata, id_var, responses)
  }
  object$X_pi <- model_matrix_initial_formula(
    initial_formula, newdata, object$n_sequences,
    object$n_states, id_var,
    X_mean = attr(object$X_pi, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_pi, "R_inv")
  )$X
  object$X_A <- model_matrix_transition_formula(
    transition_formula, newdata, object$n_sequences, object$length_of_sequences,
    object$n_states, id_var, time_var, object$sequence_lengths,
    X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_A, "R_inv")
  )$X
  object$X_B <- model_matrix_emission_formula(
    emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    id_var, time_var, object$sequence_lengths,
    X_mean = attr(object$X_B, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_B, "R_inv"), 
    autoregression = !is.null(object$autoregression_formula)
  )$X
  object$data <- newdata
  object
}
#' @rdname update_nhmm
#' @export
update.mnhmm <- function(object, newdata, ...) {
  # avoid CRAN check warning due to NSE
  tmax <- y <- NULL
  cluster_formula <- object$cluster_formula
  attr(cluster_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(cluster_formula), object$data
  )
  initial_formula <- object$initial_formula
  attr(initial_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(initial_formula), object$data
  )
  transition_formula <- object$transition_formula
  attr(transition_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(transition_formula), object$data
  )
  emission_formula <- object$emission_formula
  attr(emission_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(emission_formula), object$data
  )
  
  time_var <- object$time_variable
  id_var <- object$id_variable
  responses <- object$responses
  newdata <- .check_data(newdata, id_var, time_var, responses)
  ids <- unique(newdata[[id_var]])
  object$n_sequences <- length(ids)
  times <- unique(newdata[[time_var]])
  object$length_of_sequences <- length(times)
  n_obs <- sum(!is.na(newdata[, y, env = list(y = I(responses))])) / object$n_channels
  attr(object, "nobs") <- n_obs
  newdata[, tmax := max(time_var), by = id_var, 
          env = list(id_var = id_var, time_var = time_var)]
  newdata <- fill_time(newdata, id_var, time_var)
  object$sequence_lengths <- newdata[, 
                                     sum(time_var <= tmax, na.rm = TRUE), 
                                     by = id_var, 
                                     env = list(
                                       id_var = id_var, 
                                       time_var = time_var
                                     )]$V1
  newdata[, tmax := NULL]
  if (inherits(object, "fanhmm") && !is.null(object$autoregression_formula)) {
    newdata[[paste0("lag_", responses)]] <- 
      group_lag(newdata, id_var, responses)
  }
  object$X_omega <- model_matrix_cluster_formula(
    cluster_formula, newdata, object$n_sequences, object$n_clusters,
    id_var, X_mean = attr(object$X_omega, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_omega, "R_inv")
  )$X
  object$X_pi <- model_matrix_initial_formula(
    initial_formula, newdata, object$n_sequences, object$n_states, 
    id_var, X_mean = attr(object$X_pi, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_pi, "R_inv")
  )$X
  object$X_A <- model_matrix_transition_formula(
    transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, id_var, time_var, 
    object$sequence_lengths,
    X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_A, "R_inv")
  )$X
  object$X_B <- model_matrix_emission_formula(
    emission_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, object$n_symbols, 
    id_var, time_var, object$sequence_lengths,
    X_mean = attr(object$X_B, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_B, "R_inv"), 
    autoregression = !is.null(object$autoregression_formula)
  )$X
  object$data <- newdata
  object
}

update_W_for_fanhmm <- function(object) {
  transition_formula <- object$transition_formula
  attr(transition_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(transition_formula), object$data
  )
  emission_formula <- object$emission_formula
  attr(emission_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(emission_formula), object$data
  )
  responses <- object$responses
  newdata <- object$data
  W_A <- W_B <- vector("list", object$n_symbols)
  for (i in seq_len(object$n_symbols)) {
    newdata[[responses]][] <- object$symbol_names[i]
    newdata[[paste0("lag_", responses)]][] <- object$symbol_names[i]
    
    W_A[[i]] <- model_matrix_transition_formula(
      transition_formula, newdata, object$n_sequences, 
      object$length_of_sequences, object$n_states, object$id_variable, 
      object$time_variable, object$sequence_lengths,
      X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
      scale = TRUE, R_inv = attr(object$X_A, "R_inv")
    )$X
    W_B[[i]] <- model_matrix_emission_formula(
      emission_formula, newdata, object$n_sequences, 
      object$length_of_sequences, object$n_states, object$n_symbols, 
      object$id_variable, object$time_variable, 
      object$sequence_lengths,
      X_mean = attr(object$X_B, "X_mean"), check = FALSE, 
      scale = TRUE, R_inv = attr(object$X_B, "R_inv"), 
      autoregression = !is.null(object$autoregression_formula)
    )$X
  }
  list(W_A = W_A, W_B = W_B)
}
