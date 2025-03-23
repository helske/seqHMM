#' Does covariate values vary per ID?
#' @noRd
iv_X <- function(X) {
  all_same <- TRUE
  for (i in seq_len(dim(X)[3])) {
    all_same <- all_same && all(X[, , i] == X[, , 1])
    if (!all_same) break
  }
  !all_same #dim(unique(X, MARGIN = 3))[3] > 1L
}
#' Does covariate values vary in time?
#' @noRd
tv_X <- function(X) {
  all_same <- TRUE
  for (i in seq_len(dim(X)[2])) {
    all_same <- all_same && all(X[, i, ] == X[, 1, ])
    if (!all_same) break
  }
  !all_same #dim(unique(X, MARGIN = 2))[2] > 1L
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
model_matrix_cluster_formula <- function(formula, data, n_sequences, n_clusters, 
                                         id_var, scale = TRUE, X_mean = TRUE, 
                                         R_inv = NULL, check = TRUE) {
  icpt_only <- intercept_only(formula)
  if (icpt_only) {
    n_pars <- n_clusters - 1L
    X <- matrix(1, n_sequences, 1)
    coef_names <- "(Intercept)"
    iv <- FALSE
    X_mean <- NULL
    R_inv <- NULL
  } else {
    data <- data[, .SD[1], by = id_var]
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    missing_values <- which(!stats::complete.cases(X))
    stopifnot_(
      length(missing_values) == 0L,
      c(
        "Missing cases are not allowed in covariates of `cluster_formula`.",
        "Use {.fn stats::complete.cases} to detect them, then fix or impute them.",
        paste0(
          "First missing value found for ID ",
          "{data[missing_values[1], c(id_var)]}."
        )
      )
    )
    n_pars <- (n_clusters - 1L) * ncol(X)
    if (check) .check_identifiability(X, "cluster")
    
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    
    if (scale) {
      X_scaled <- scale(X[, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        scale_n <- sqrt(nrow(X_scaled) - 1)
        R_inv <- solve(qr.R(qr_out)[, qr_out$pivot] / scale_n)
        X[, cols] <- qr.Q(qr_out) * scale_n
      } else {
        X[, cols] <- X_scaled %*% R_inv #t(solve(t(R), t(X_scaled)))
      }
    } else {
      X_mean <- rep(0, length(cols))
      R_inv <- diag(length(cols))
    }
  }
  X <- t(X)
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "icpt_only") <- icpt_only
  
  list(formula = formula, n_pars = n_pars, X = X)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_initial_formula <- function(formula, data, n_sequences, n_states, 
                                         id_var, scale = TRUE, X_mean = TRUE, 
                                         R_inv = NULL, check = TRUE) {
  icpt_only <- intercept_only(formula)
  if (icpt_only) {
    n_pars <- n_states - 1L
    X <- matrix(1, n_sequences, 1)
    coef_names <- "(Intercept)"
    X_mean <- NULL
    R_inv <- NULL
  } else {
    data <- data[, .SD[1], by = id_var]
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    missing_values <- which(!stats::complete.cases(X))
    stopifnot_(
      length(missing_values) == 0L,
      c(
        "Missing cases are not allowed in covariates of `initial_formula`.",
        "Use {.fn stats::complete.cases} to detect them, then fix or impute them.",
        paste0(
          "First missing value found for ID ",
          "{data[missing_values[1], c(id_var)]}."
        )
      )
    )
    n_pars <- (n_states - 1L) * ncol(X)
    if (check) .check_identifiability(X, "initial")
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    if (scale) {
      X_scaled <- scale(X[, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        scale_n <- sqrt(nrow(X_scaled) - 1)
        R_inv <- solve(qr.R(qr_out)[, qr_out$pivot] / scale_n)
        X[, cols] <- qr.Q(qr_out) * scale_n
      } else {
        X[, cols] <- X_scaled %*% R_inv #t(solve(t(R), t(X_scaled)))
      }
    } else {
      X_mean <- rep(0, length(cols))
      R_inv <- diag(length(cols))
    }
  }
  X <- t(X)
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "icpt_only") <- icpt_only
  
  list(formula = formula, n_pars = n_pars, X = X)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_transition_formula <- function(formula, data, n_sequences,
                                            length_of_sequences, n_states,
                                            id_var, time_var, sequence_lengths, 
                                            scale = TRUE, X_mean = TRUE, 
                                            R_inv = NULL, check = TRUE) {
  # avoid CRAN check warning due to NSE
  time <- NULL
  icpt_only <- intercept_only(formula)
  if (icpt_only) {
    n_pars <-  n_states * (n_states - 1L)
    X <- array(1, c(1L, length_of_sequences, n_sequences))
    coef_names <- "(Intercept)"
    iv <- tv <- FALSE
    missing_values <- integer(0)
    X_mean <- NULL
    R_inv <- NULL
  } else {
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    complete <- stats::complete.cases(X)
    missing_values <- which(!complete)
    # A_1 is not used anywhere, so we can allow missing values in the first time point
    missing_values <- setdiff(missing_values, which(data[[time_var]] == min(data[[time_var]])))
    if (length(missing_values) > 0 && check) {
      idx <- sequence_lengths[match(data[[id_var]], unique(data[[id_var]]))]
      ends <- data[idx, time, env = list(time = time_var)]
      stopifnot_(
        all(z <- data[missing_values, time, env = list(time = time_var)] > ends[missing_values]),
        c(
          paste0(
            "Missing cases are not allowed in covariates of ",
            "{.arg transition_formula}, unless they correspond to void ",
            "response values at the end of the sequences.",
            "Use {.fn stats::complete.cases} to detect them, then fix or impute them."
          ),
          paste0(
            "First missing value found for ID ",
            "{data[missing_values, c(id_var)][which(!z)[1]]} at time point ",
            "{data[missing_values, c(time_var)][which(!z)[1]]}."
          )
        )
      )
    }
    if (check) .check_identifiability(X[complete, ], "transition")
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    if (scale) {
      X_scaled <- scale(X[complete, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        scale_n <- sqrt(nrow(X_scaled) - 1)
        R_inv <- solve(qr.R(qr_out)[, qr_out$pivot] / scale_n)
        X[complete, cols] <- qr.Q(qr_out) * scale_n
      } else {
        X[complete, cols] <- X_scaled %*% R_inv #t(solve(t(R), t(X_scaled)))
      }
    } else {
      X_mean <- rep(0, length(cols))
      R_inv <- diag(length(cols))
    }
    X <- t(X)
    dim(X) <- c(nrow(X), length_of_sequences, n_sequences)
    n_pars <- n_states * (n_states - 1L) * nrow(X)
    missing_values <- which(is.na(X))
    # Replace NAs in void cases and t = 1 with zero
    X[is.na(X)] <- 0
    iv <- iv_X(X[, -1L, , drop = FALSE])
    tv <- tv_X(X[, -1L, , drop = FALSE])
  }
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  attr(X, "tv") <- tv
  attr(X, "icpt_only") <- icpt_only
  attr(X, "missing") <- missing_values
  
  list(formula = formula, n_pars = n_pars, X = X)
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_emission_formula <- function(formula, data, n_sequences, 
                                          length_of_sequences, n_states,
                                          n_symbols,
                                          id_var, time_var, sequence_lengths, 
                                          scale = TRUE, X_mean = TRUE, 
                                          R_inv = NULL, check = TRUE,
                                          autoregression = FALSE) {
  # avoid CRAN check warning due to NSE
  time <- NULL
  icpt_only <- intercept_only(formula)
  if (icpt_only) {
    n_pars <-  sum(n_states * (n_symbols - 1L))
    X <- array(1, c(1L, length_of_sequences, n_sequences))
    coef_names <- "(Intercept)"
    iv <- tv <- FALSE
    missing_values <- integer(0)
    X_mean <- NULL
    R_inv <- NULL
  } else {
    X <- stats::model.matrix.lm(
      formula, 
      data = data, 
      na.action = stats::na.pass
    )
    
    complete <- stats::complete.cases(X)
    missing_values <- which(!complete)
    if (autoregression) {
      # first observation is fixed so missing (lagged) values do not matter
      missing_values <- setdiff(missing_values, which(data[[time_var]] == min(data[[time_var]])))
    } 
    if (length(missing_values) > 0 && check) {
      idx <- sequence_lengths[match(data[[id_var]], unique(data[[id_var]]))]
      ends <- data[idx, time, env = list(time = time_var)]
      stopifnot_(
        all(z <- data[missing_values, time, env = list(time = time_var)] > ends[missing_values]),
        c(paste0(
          "Missing cases are not allowed in covariates of ",
          "{.arg emission_formula}, unless they correspond to missing ",
          "responses at the end of the sequences. ",
          "Use {.fn stats::complete.cases} to detect them, then fix or impute them. "
        ),
        paste0(
          "First missing value found for ID ",
          "{data[missing_values, c(id_var)][which(!z)[1]]} at time point ",
          "{data[missing_values, c(time_var)][which(!z)[1]]}."
        )
        )
      )
    }
    if (check) .check_identifiability(X[complete, ], "emission")
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    if (scale) {
      X_scaled <- scale(X[complete, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        scale_n <- sqrt(nrow(X_scaled) - 1)
        R_inv <- solve(qr.R(qr_out)[, qr_out$pivot] / scale_n)
        X[complete, cols] <- qr.Q(qr_out) * scale_n
      } else {
        X[complete, cols] <- X_scaled %*% R_inv #t(solve(t(R), t(X_scaled)))
      }
    } else {
      X_mean <- rep(0, length(cols))
      R_inv <- diag(length(cols))
    }
    X <- t(X)
    dim(X) <- c(nrow(X), length_of_sequences, n_sequences)
    n_pars <- sum(n_states * (n_symbols - 1L) * nrow(X))
    missing_values <- which(is.na(X))
    # Replace NAs in void cases with zero
    X[is.na(X)] <- 0
    iv <- iv_X(X)
    tv <- tv_X(X)
  }
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  attr(X, "tv") <- tv
  attr(X, "icpt_only") <- icpt_only
  attr(X, "missing") <- missing_values
  
  list(formula = formula, n_pars = n_pars, X = X)
}
