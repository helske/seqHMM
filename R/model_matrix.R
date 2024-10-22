#' Does covariate values vary per ID?
#' @noRd
iv_X <- function(X) {
  dim(unique(X, MARGIN = 3))[3] > 1L
}
#' Does covariate values vary in time?
#' @noRd
tv_X <- function(X) {
  dim(unique(X, MARGIN = 2))[2] > 1L
}

# Function to check uniqueness along the N dimension
check_unique_N <- function(arr) {
  # Flatten along the T and K dimensions
  flattened_N <- apply(arr, 2, function(x) as.vector(x))
  # Check for unique rows
  length(unique(as.data.frame(t(flattened_N)))) > 1
}

#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_initial_formula <- function(formula, data, n_sequences, 
                                         length_of_sequences, n_states, 
                                         time, id, X_mean, X_sd) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    n_pars <- n_states - 1L
    X <- matrix(1, n_sequences, 1)
    coef_names <- "(Intercept)"
    iv <- FALSE
    X_mean <- NULL
    X_sd <- NULL
  } else {
    first_time_point <- min(data[[time]])
    data <- data[data[[time]] == first_time_point, ]
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    cols <- which(colnames(X) != "(Intercept)")
    if (missing(X_mean)) {
      X_mean <- X_sd <- TRUE
    }
    X_scaled <- scale(X[, cols], X_mean, X_sd)
    X[, cols] <- X_scaled
    X_mean <- attr(X_scaled, "scaled:center")
    X_sd <- attr(X_scaled, "scaled:scale")
    
    missing_values <- which(!complete.cases(X))
    stopifnot_(
      length(missing_values) == 0L,
      c(
        "Missing cases are not allowed in covariates of `initial_formula`.",
        "Use {.fn complete.cases} to detect them, then fix or impute them.",
        paste0(
          "First missing value found for ID ",
          "{data[missing_values[1], id]}."
        )
      )
    )
    iv <- nrow(unique(X)) > 1L
    coef_names <- colnames(X)
    n_pars <- (n_states - 1L) * ncol(X)
  }
  X <- t(X)
  attr(X, "X_mean") <- X_mean
  attr(X, "X_sd") <- X_sd
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  list(formula = formula, n_pars = n_pars, X = X)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_transition_formula <- function(formula, data, n_sequences,
                                            length_of_sequences, n_states,
                                            time, id, sequence_lengths, 
                                            X_mean, X_sd) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    n_pars <-  n_states * (n_states - 1L)
    X <- array(1, c(1L, length_of_sequences, n_sequences))
    coef_names <- "(Intercept)"
    iv <- tv <- FALSE
    missing_values <- integer(0)
    X_mean <- NULL
    X_sd <- NULL
  } else {
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    cols <- which(colnames(X) != "(Intercept)")
    if (missing(X_mean)) {
      X_mean <- X_sd <- TRUE
    }
    X_scaled <- scale(X[, cols], X_mean, X_sd)
    X[, cols] <- X_scaled
    X_mean <- attr(X_scaled, "scaled:center")
    X_sd <- attr(X_scaled, "scaled:scale")
    
    missing_values <- which(!complete.cases(X))
    if (length(missing_values) > 0) {
      ends <- sequence_lengths[match(data[[id]], unique(data[[id]]))]
      stopifnot_(
        all(z <- data[missing_values, time] > ends[missing_values]),
        c(
          paste0(
            "Missing cases are not allowed in covariates of ",
            "{.arg transition_formula}, unless they correspond to void ",
            "response values at the end of the sequences.",
            "Use {.fn complete.cases} to detect them, then fix or impute them."
          ),
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
    n_pars <- n_states * (n_states - 1L) * dim(X)[3]
    iv <- iv_X(X)
    tv <- tv_X(X)
    X <- aperm(X, c(3, 1, 2))
    missing_values <- which(is.na(X))
    # Replace NAs in void cases with zero
    X[is.na(X)] <- 0
  }
  attr(X, "X_mean") <- X_mean
  attr(X, "X_sd") <- X_sd
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  attr(X, "tv") <- tv
  attr(X, "missing") <- missing_values
  list(formula = formula, n_pars = n_pars, X = X)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_emission_formula <- function(formula, data, n_sequences, 
                                          length_of_sequences, n_states,
                                          n_symbols, n_channels,
                                          time, id, sequence_lengths, 
                                          X_mean, X_sd) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    n_pars <-  sum(n_states * (n_symbols - 1L))
    X <- array(1, c(1L, length_of_sequences, n_sequences))
    coef_names <- "(Intercept)"
    iv <- tv <- FALSE
    missing_values <- integer(0)
    X_mean <- NULL
    X_sd <- NULL
  } else {
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    cols <- which(colnames(X) != "(Intercept)")
    if (missing(X_mean)) {
      X_mean <- X_sd <- TRUE
    }
    X_scaled <- scale(X[, cols], X_mean, X_sd)
    X[, cols] <- X_scaled
    X_mean <- attr(X_scaled, "scaled:center")
    X_sd <- attr(X_scaled, "scaled:scale")
    
    missing_values <- which(!complete.cases(X))
    if (length(missing_values) > 0) {
      ends <- sequence_lengths[match(data[[id]], unique(data[[id]]))]
      stopifnot_(
        all(z <- data[missing_values, time] > ends[missing_values]),
        c(paste0(
          "Missing cases are not allowed in covariates of ",
          "{.arg emission_formula}, unless they correspond to missing ",
          "void reponses at the end of the sequences. ",
          "Use {.fn complete.cases} to detect them, then fix or impute them. ",
          "Note that the missing covariates in {.arg emission_formula} ",
          "corresponding to time points where response variables are also ",
          "missing can be set to arbitrary value."
        ),
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
    n_pars <- sum(n_states * (n_symbols - 1L) * dim(X)[3])
    iv <- iv_X(X)
    tv <- tv_X(X)
    X <- aperm(X, c(3, 1, 2))
    missing_values <- which(is.na(X))
    # Replace NAs in void cases with zero
    X[is.na(X)] <- 0
  }
  attr(X, "X_mean") <- X_mean
  attr(X, "X_sd") <- X_sd
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  attr(X, "tv") <- tv
  attr(X, "missing") <- missing_values
  list(formula = formula, n_pars = n_pars, X = X)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_cluster_formula <- function(formula, data, n_sequences, n_clusters, 
                                         time, id, X_mean, X_sd) {
  icp_only <- intercept_only(formula)
  if (icp_only) {
    n_pars <- n_clusters - 1L
    X <- matrix(1, n_sequences, 1)
    coef_names <- "(Intercept)"
    iv <- FALSE
    X_mean <- NULL
    X_sd <- NULL
  } else {
    first_time_point <- min(data[[time]])
    X <- stats::model.matrix.lm(
      formula, 
      data = data[data[[time]] == first_time_point, ], 
      na.action = stats::na.pass
    )
    cols <- which(colnames(X) != "(Intercept)") #always first column(?)
    if (missing(X_mean)) {
      X_mean <- X_sd <- TRUE
    }
    X_scaled <- scale(X[, cols], X_mean, X_sd)
    X[, cols] <- X_scaled
    X_mean <- attr(X_scaled, "scaled:center")
    X_sd <- attr(X_scaled, "scaled:scale")
    
    missing_values <- which(!complete.cases(X))
    stopifnot_(
      length(missing_values) == 0L,
      c(
        "Missing cases are not allowed in covariates of `cluster_formula`.",
        "Use {.fn complete.cases} to detect them, then fix or impute them.",
        paste0(
          "First missing value found for ID ",
          "{data[missing_values[1], id]}."
        )
      )
    )
    iv <- nrow(unique(X)) > 1L
    coef_names <- colnames(X)
    n_pars <- (n_clusters - 1L) * ncol(X)
  }
  X <- t(X)
  attr(X, "X_mean") <- X_mean
  attr(X, "X_sd") <- X_sd
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  list(formula = formula, n_pars = n_pars, X = X)
}
