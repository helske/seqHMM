#' Update Covariate Values of NHMM
#' 
#' This function can be used to replace original covariate values of NHMMs. 
#' The responses, model formulae and estimated coefficients are not altered.
#' @param object An object of class `nhmm` or `mnhmm`.
#' @param newdata A data frame containing the new covariate values.
#' @param drop_levels if `TRUE` (default), drops unused factor levels from 
#' `newdata` before creating design matrices.
#' @param ... Ignored.
#' @rdname update_nhmm
#' @export
update.nhmm <- function(object, newdata, drop_levels = TRUE, ...) {
  # avoid CRAN check warning due to NSE
  .Ti <- y <- id <- NULL
  initial_formula <- object$initial_formula
  attr(initial_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(initial_formula), object$data
  )
  transition_formula <- object$transition_formula
  attr(transition_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(transition_formula), object$data
  )
  emission_formula <- object$emission_formula
  for (i in seq_len(object$n_channels)) {
    attr(emission_formula[[i]], "xlevels") <- stats::.getXlevels(
      stats::terms(emission_formula[[i]]), object$data
    )
  }
  time_var <- object$time_variable
  id_var <- object$id_variable
  responses <- object$responses
  newdata <- .check_data(newdata, id_var, time_var, responses)
  newdata <- fill_time(newdata, id_var, time_var)
  object$sequence_lengths <- newdata[, .Ti[1], by = id_var, 
                                     env = list(id_var = id_var), 
                                     showProgress = FALSE]$V1
  newdata[, .Ti := NULL]
  
  if (inherits(object, "fanhmm")) {
    for (y in responses) {
      lag_obs <- paste0("lag_", y)
      y1 <- newdata[[y]][1] #the value is not used anywhere
      newdata[, lag_obs := shift(y, type = "lag", fill = y1), by = id, 
              env = list(id = id_var, y = y, lag_obs = lag_obs, y1 = y1), 
              showProgress = FALSE]
    }
    if (length(object$autoregression) > 0L && identical(object$prior_obs, 0L)) {
      .idx <- setdiff(
        seq_row(newdata), 
        cumsum(c(1, head(object$sequence_lengths, -1)))
      )
      newdata <- newdata[.idx]
      object$sequence_lengths <- object$sequence_lengths - 1L
      msg <- paste0(
        "The model contains lagged responses in emission formula and first time", 
        " point was fixed ({.arg prior_obs} = {.val fixed}). ",
        " {.arg newdata} should contain also this time point."
      )
    }
  } else {
    msg <- NULL
  }
  if (drop_levels) {
    setdroplevels(newdata)
  }
  times_old <- unique(
    object$data[[time_var]], 
    nmax = object$length_of_sequences
  )
  times_new <- unique(newdata[[time_var]])
  stopifnot_(
    min(times_new) == min(times_old),
    c("Earliest time point in {.arg newdata} should match the earliest time point 
    in the original data.",
      i = msg)
  )
  stopifnot_(
    min(diff(times_new)) == min(diff(times_old)),
    "Time resolution in {.arg newdata} should match the resolution in the 
    original data."
  )
  
  object$length_of_sequences <- n_unique(newdata[[time_var]])
  object$n_sequences <- n_unique(newdata[[id_var]])
  n_obs <- sum(!is.na(newdata[, y, env = list(y = I(responses))])) / object$n_channels
  attr(object, "nobs") <- n_obs
  
  object$X_pi <- model_matrix_initial_formula(
    initial_formula, newdata, object$n_sequences,
    object$n_states, id_var,
    X_mean = attr(object$X_pi, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_pi, "R_inv")
  )
  object$X_A <- model_matrix_transition_formula(
    transition_formula, newdata, object$n_sequences, object$length_of_sequences,
    object$n_states, id_var, time_var, object$sequence_lengths,
    X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_A, "R_inv")
  )
  object$X_B <- stats::setNames(
    lapply(responses, \(y) {
      model_matrix_emission_formula(
        emission_formula[[y]], newdata, object$n_sequences, 
        object$length_of_sequences, object$n_states, 
        object$n_symbols[y], id_var, time_var, object$sequence_lengths, 
        X_mean = attr(object$X_B[[y]], "X_mean"), check = FALSE, 
        scale = TRUE, R_inv = attr(object$X_B[[y]], "R_inv")
      )
    }),
    responses
  )
  if (inherits(object, "fanhmm") && !identical(object$prior_obs, 0L)) {
    object$W_X_B <- create_W_X_B(
      newdata, id_var, time_var, object$symbol_names, object$n_sequences, 
      emission_formula, object$n_states, object$X_B
    )
  }
  object$data <- newdata
  object
}
#' @rdname update_nhmm
#' @export
update.mnhmm <- function(object, newdata, drop_levels = TRUE, ...) {
  # avoid CRAN check warning due to NSE
  .Ti <- y <- id <- NULL
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
  for (i in seq_len(object$n_channels)) {
    attr(emission_formula[[i]], "xlevels") <- stats::.getXlevels(
      stats::terms(emission_formula[[i]]), object$data
    )
  }
  time_var <- object$time_variable
  id_var <- object$id_variable
  responses <- object$responses
  newdata <- .check_data(newdata, id_var, time_var, responses)
  newdata <- fill_time(newdata, id_var, time_var)
  object$sequence_lengths <- newdata[, .Ti[1], by = id_var, 
                                     env = list(id_var = id_var), 
                                     showProgress = FALSE]$V1
  newdata[, .Ti := NULL]
  if (inherits(object, "fanhmm")) {
    for (y in responses) {
      lag_obs <- paste0("lag_", y)
      y1 <- newdata[[y]][1] #the value is not used anywhere
      newdata[, lag_obs := shift(y, type = "lag", fill = y1), by = id, 
              env = list(id = id_var, y = y, lag_obs = lag_obs, y1 = y1), 
              showProgress = FALSE]
    }
    if (length(object$autoregression) > 0L && identical(object$prior_obs, 0L)) {
      .idx <- setdiff(
        seq_row(newdata), 
        cumsum(c(1, head(object$sequence_lengths, -1)))
      )
      newdata <- newdata[.idx]
      object$sequence_lengths <- object$sequence_lengths - 1L
      msg <- paste0(
        "The model contains lagged responses in emission formula and first time", 
        " point was fixed ({.arg prior_obs} = {.val fixed}). ",
        " {.arg newdata} should contain also this time point."
      )
    }
  } else {
    msg <- NULL
  }
  if (drop_levels) {
    setdroplevels(newdata)
  }
  times_old <- unique(
    object$data[[time_var]], 
    nmax = object$length_of_sequences
  )
  times_new <- unique(newdata[[time_var]])
  stopifnot_(
    min(times_new) == min(times_old),
    c("Earliest time point in {.arg newdata} should match the earliest time point 
    in the original data.",
      i = msg)
  )
  stopifnot_(
    min(diff(times_new)) == min(diff(times_old)),
    "Time resolution in {.arg newdata} should match the resolution in the 
    original data."
  )
  object$length_of_sequences <- n_unique(newdata[[time_var]])
  object$n_sequences <- n_unique(newdata[[id_var]])
  n_obs <- sum(!is.na(newdata[, y, env = list(y = I(responses))])) / object$n_channels
  attr(object, "nobs") <- n_obs
  
  object$X_omega <- model_matrix_cluster_formula(
    cluster_formula, newdata, object$n_sequences, object$n_clusters,
    id_var, X_mean = attr(object$X_omega, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_omega, "R_inv")
  )
  object$X_pi <- model_matrix_initial_formula(
    initial_formula, newdata, object$n_sequences, object$n_states, 
    id_var, X_mean = attr(object$X_pi, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_pi, "R_inv")
  )
  object$X_A <- model_matrix_transition_formula(
    transition_formula, newdata, object$n_sequences, 
    object$length_of_sequences, object$n_states, id_var, time_var, 
    object$sequence_lengths,
    X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
    scale = TRUE, R_inv = attr(object$X_A, "R_inv")
  )
  object$X_B <- stats::setNames(
    lapply(responses, \(y) {
      model_matrix_emission_formula(
        emission_formula[[y]], newdata, object$n_sequences, 
        object$length_of_sequences, object$n_states, 
        object$n_symbols[y], id_var, time_var, object$sequence_lengths, 
        X_mean = attr(object$X_B[[y]], "X_mean"), check = FALSE, 
        scale = TRUE, R_inv = attr(object$X_B[[y]], "R_inv")
      )
    }),
    responses
  )
  if (length(object$autoregression) > 0L && !identical(object$prior_obs, 0L)) {
    object$W_X_B <- create_W_X_B(
      newdata, id_var, time_var, object$symbol_names, object$n_sequences, 
      emission_formula, object$n_states, object$X_B
    )
  }
  object$data <- newdata
  object
}

update_W_for_fanhmm <- function(object) {
  
  transition_formula <- object$transition_formula
  attr(transition_formula, "xlevels") <- stats::.getXlevels(
    stats::terms(transition_formula), object$data
  )
  emission_formula <- object$emission_formula
  for (i in seq_len(object$n_channels)) {
    attr(emission_formula[[i]], "xlevels") <- stats::.getXlevels(
      stats::terms(emission_formula[[i]]), object$data
    )
  }
  responses <- object$responses
  newdata <- object$data
  C <- length(responses)
  y <- expand.grid(object$symbol_names)
  W_A <- W_B <- vector("list", nrow(y))
  for (i in seq_along(W_A)) {
    for (j in seq_len(C)) {
      set(newdata, j = paste0("lag_", responses[j]), value = y[i, j])
    }
    W_A[[i]] <- model_matrix_transition_formula(
      transition_formula, newdata, object$n_sequences, 
      object$length_of_sequences, object$n_states, object$id_variable, 
      object$time_variable, object$sequence_lengths,
      X_mean = attr(object$X_A, "X_mean"), check = FALSE, 
      scale = TRUE, R_inv = attr(object$X_A, "R_inv")
    )
    W_B[[i]] <- vector("list", C)
    for (j in seq_len(C)) {
      W_B[[i]][[j]] <- model_matrix_emission_formula(
        emission_formula[[j]], newdata, object$n_sequences, 
        object$length_of_sequences, object$n_states, object$n_symbols[[j]], 
        object$id_variable, object$time_variable, 
        object$sequence_lengths,
        X_mean = attr(object$X_B[[j]], "X_mean"), check = FALSE, 
        scale = TRUE, R_inv = attr(object$X_B[[j]], "R_inv")
      )
    }
  }
  list(W_A = W_A, W_B = W_B)
}
