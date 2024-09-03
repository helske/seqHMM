#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_initial_formula <- function(formula, data, n_sequences, 
                                         length_of_sequences, n_states, 
                                         time, id) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    type <- "c"
    n_pars <- n_states - 1L
    X <- matrix(1, n_sequences, 1)
    coef_names <- "(Intercept)"
  } else {
    first_time_point <- min(data[[time]])
    X <- stats::model.matrix.lm(
      formula, 
      data = data[data[[time]] == first_time_point, ], 
      na.action = stats::na.pass
    )
    missing_values <- which(!complete.cases(X))
    stopifnot_(
      length(missing_values) == 0,
      c(
        "Missing cases are not allowed in covariates of `initial_formula`.",
        "Use {.fn complete.cases} to detect them, then fix or impute them.",
        paste0(
          "First missing value found for ID ",
          "{data[missing_values[1], id]}."
        )
      )
    )
    coef_names <- colnames(X)
    type <- "v"
    n_pars <- (n_states - 1L) * ncol(X)
  }
  list(formula = formula, type = type, n_pars = n_pars, X = X, 
       coef_names = coef_names)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_transition_formula <- function(formula, data, n_sequences,
                                            length_of_sequences, n_states,
                                            time, id, sequence_lengths) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    type <- "c"
    n_pars <-  n_states * (n_states - 1L)
    X <- array(1, c(length_of_sequences, n_sequences, 1L))
    coef_names <- "(Intercept)"
  } else {
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    missing_values <- which(!complete.cases(X))
    if (length(missing_values) > 0) {
      ends <- sequence_lengths[match(data[[id]], unique(data[[id]]))]
      stopifnot_(
        all(z <- data[missing_values, time] <= ends[missing_values]),
        c(
          "Missing cases are not allowed in covariates of `transition_formula`.",
          "Use {.fn complete.cases} to detect them, then fix or impute them.",
          paste0(
            "First missing value found for ID ",
            "{data[missing_values, id][which(!z)[1]]} at time point ",
            "{data[missing_values, time][which(!z)[1]]}."
          )
        )
      )
    }
    coef_names <- colnames(X)
    dim(X) <- c(length_of_sequences, n_sequences, ncol(X))
    type <- "v"
    n_pars <- n_states * (n_states - 1L) * dim(X)[3]
  }
  list(formula = formula, type = type, n_pars = n_pars, X = X, 
       coef_names = coef_names)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_emission_formula <- function(formula, data, n_sequences, 
                                          length_of_sequences, n_states,
                                          n_symbols, n_channels,
                                          time, id, sequence_lengths) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    type <- "c"
    n_pars <-  n_channels * n_states * (n_symbols - 1L)
    X <- array(1, c(length_of_sequences, n_sequences, 1L))
    coef_names <- "(Intercept)"
  } else {
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    missing_values <- which(!complete.cases(X))
    if (length(missing_values) > 0) {
      ends <- sequence_lengths[match(data[[id]], unique(data[[id]]))]
      stopifnot_(
        all(z <- data[missing_values, time] <= ends[missing_values]),
        c(
          "Missing cases are not allowed in covariates of `emission_formula`.",
          "Use {.fn complete.cases} to detect them, then fix or impute them.",
          paste0(
            "First missing value found for ID ",
            "{data[missing_values, id][which(!z)[1]]} at time point ",
            "{data[missing_values, time][which(!z)[1]]}."
          )
        )
      )
    }
    coef_names <- colnames(X)
    dim(X) <- c(length_of_sequences, n_sequences, ncol(X))
    type <- "v"
    n_pars <- n_channels * n_states * (n_symbols - 1L) * dim(X)[3]
  }
  list(formula = formula, type = type, n_pars = n_pars, X = X, 
       coef_names = coef_names)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_cluster_formula <- function(formula, data, n_sequences, n_clusters, 
                                         time, id) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    type <- "c"
    n_pars <- n_clusters - 1L
    X <- matrix(1, n_sequences, 1)
    coef_names <- "(Intercept)"
  } else {
    first_time_point <- min(data[[time]])
    X <- stats::model.matrix.lm(
      formula, 
      data = data[data[[time]] == first_time_point, ], 
      na.action = stats::na.pass
    )
    missing_values <- which(!complete.cases(X))
    stopifnot_(
      length(missing_values) == 0,
      c(
        "Missing cases are not allowed in covariates of `cluster_formula`.",
        "Use {.fn complete.cases} to detect them, then fix or impute them.",
        paste0(
          "First missing value found for ID ",
          "{data[missing_values[1], id]}."
        )
      )
    )
    coef_names <- colnames(X)
    type <- "v"
    n_pars <- (n_clusters - 1L) * ncol(X)
  }
  list(formula = formula, type = type, n_pars = n_pars, X = X, 
       coef_names = coef_names)
}
