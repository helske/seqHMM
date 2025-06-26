#' Does covariate values vary in time?
#' @noRd
tv_X <- function(X, cols) {
  any(fndistinct(X[, cols, env = list(cols = I(cols))], g = X$id) > 1)
}
#' Does covariate values vary per ID?
#' @noRd
iv_X <- function(X, cols) {
  any(fndistinct(X[, cols, env = list(cols = I(cols))], g = X$time) > 1)
}

#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_cluster_formula <- function(formula, data, n_sequences, n_clusters, 
                                         id_var, scale = TRUE, X_mean = TRUE, 
                                         R_inv = NULL, check = TRUE) {
  icpt_only <- intercept_only(formula)
  if (icpt_only) {
    X <- matrix(1, 1, n_sequences)
    coef_names <- "(Intercept)"
    iv <- FALSE
    X_mean <- NULL
    R_inv <- NULL
  } else {
    data <- data[, .SD[1], by = id_var, showProgress = FALSE]
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
          "{data[missing_values[1], id, env = list(id = id_var)]}."
        )
      )
    )
    if (check && !scale && is.null(R_inv)) .check_identifiability(X, "cluster")
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    if (scale) {
      X_scaled <- scale(X[, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        if (check) .check_identifiability(X_scaled, "cluster", qr_out)
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
    X <- t(X)
  }
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "icpt_only") <- icpt_only
  X
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_initial_formula <- function(formula, data, n_sequences, n_states, 
                                         id_var, scale = TRUE, X_mean = TRUE, 
                                         R_inv = NULL, check = TRUE) {
  icpt_only <- intercept_only(formula)
  if (icpt_only) {
    X <- matrix(1, 1, n_sequences)
    coef_names <- "(Intercept)"
    X_mean <- NULL
    R_inv <- NULL
  } else {
    data <- data[, .SD[1], by = id_var, showProgress = FALSE]
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
        `x` = paste0(
          "First missing value found for ID ",
          "{data[missing_values[1], id, env = list(id = id_var)]}."
        )
      )
    )
    if (check && !scale && is.null(R_inv)) .check_identifiability(X, "initial")
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    if (scale) {
      X_scaled <- scale(X[, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        if (check) .check_identifiability(X_scaled, "initial", qr_out)
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
    X <- t(X)
  }
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "icpt_only") <- icpt_only
  X
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
    X <- lapply(seq_len(n_sequences), \(i) matrix(1, 1, sequence_lengths[i]))
    coef_names <- "(Intercept)"
    iv <- n_unique(sequence_lengths) != 1L
    tv <- FALSE
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
      .idx <- sequence_lengths[match(data[[id_var]], unique(data[[id_var]]))]
      ends <- data[.idx, time, env = list(time = time_var)]
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
            "{data[missing_values, id, env = list(id = id_var)][which(!z)[1]]} at time point ",
            "{data[missing_values, time, env = list(time = time_var)][which(!z)[1]]}."
          )
        )
      )
    }
    if (check && !scale && is.null(R_inv)) .check_identifiability(X[complete, ], "transition")
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    if (scale) {
      X_scaled <- scale(X[complete, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        if (check) .check_identifiability(X_scaled, "transition", qr_out)
        scale_n <- sqrt(nrow(X_scaled) - 1)
        R_inv <- solve(qr.R(qr_out)[, qr_out$pivot] / scale_n)
        X[complete, cols] <- qr.Q(qr_out) * scale_n
      } else {
        X[complete, cols] <- X_scaled %*% R_inv
      }
    } else {
      X_mean <- rep(0, length(cols))
      R_inv <- diag(length(cols))
    }
    K <- ncol(X)
    X <- data.table(id = data[[id_var]], time = data[[time_var]], X, 
                       key = c("id", "time"))
    tv <- tv_X(X, coef_names[cols])
    iv <- n_unique(sequence_lengths) != 1L || iv_X(X, coef_names[cols])
    X[, time := NULL]
    X <- lapply(rsplit(X, by = ~ id), \(x) t(qM(x)))
  }
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  attr(X, "tv") <- tv
  attr(X, "icpt_only") <- icpt_only
  X
}
#' Create the Model Matrix based on NHMM Formulas
#'
#' @noRd
model_matrix_emission_formula <- function(formula, data, n_sequences, 
                                          length_of_sequences, n_states,
                                          n_symbols,
                                          id_var, time_var, sequence_lengths, 
                                          scale = TRUE, X_mean = TRUE, 
                                          R_inv = NULL, check = TRUE) {
  # avoid CRAN check warning due to NSE
  time <- NULL
  icpt_only <- intercept_only(formula)
  if (icpt_only) {
    X <- lapply(seq_len(n_sequences), \(i) matrix(1, 1, sequence_lengths[i]))
    coef_names <- "(Intercept)"
    iv <- n_unique(sequence_lengths) != 1L
    tv <- FALSE
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
    if (length(missing_values) > 0 && check) {
      .idx <- sequence_lengths[match(data[[id_var]], unique(data[[id_var]]))]
      ends <- data[.idx, time, env = list(time = time_var)]
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
          "{data[missing_values, id, env = list(id = id_var)][which(!z)[1]]} at time point ",
          "{data[missing_values, time, env = list(time = time_var)][which(!z)[1]]}."
        )
        )
      )
    }
    if (check && !scale && is.null(R_inv)) .check_identifiability(X[complete, ], "emission")
    coef_names <- colnames(X)
    cols <- which(coef_names != "(Intercept)")
    if (scale) {
      X_scaled <- scale(X[complete, cols], X_mean, FALSE)
      X_mean <- attr(X_scaled, "scaled:center")
      if (is.null(R_inv)) {
        qr_out <- qr(X_scaled)
        if (check) .check_identifiability(X_scaled, "emission", qr_out)
        scale_n <- sqrt(nrow(X_scaled) - 1)
        R_inv <- solve(qr.R(qr_out)[, qr_out$pivot] / scale_n)
        X[complete, cols] <- qr.Q(qr_out) * scale_n
      } else {
        X[complete, cols] <- X_scaled %*% R_inv
      }
    } else {
      X_mean <- rep(0, length(cols))
      R_inv <- diag(length(cols))
    }
    K <- ncol(X)
    X <- data.table(id = data[[id_var]], time = data[[time_var]], X, 
                    key = c("id", "time"))
    tv <- tv_X(X, coef_names[cols])
    iv <- n_unique(sequence_lengths) != 1L || iv_X(X, coef_names[cols])
    X[, time := NULL]
    X <- lapply(rsplit(X, by = ~ id), \(x) t(qM(x)))
  }
  attr(X, "R_inv") <- R_inv
  attr(X, "X_mean") <- X_mean
  attr(X, "coef_names") <- coef_names
  attr(X, "iv") <- iv
  attr(X, "tv") <- tv
  attr(X, "icpt_only") <- icpt_only
  X
}
#' Create the design matrix for the emissions at first time point
#' This is for marginalization over t=0
#' @noRd
create_W_X_B <- function(data, id_var, time_var, symbol_names, n_sequences,
                         emission_formula, n_states, X_B) {
  new <- old <- NULL
  combs <- expand.grid(symbol_names)
  idx <- data[, .I[1], by = id_var, env = list(id_var = id_var),
              showProgress = FALSE]$V1
  d <- data[idx]
  W_X_B <- vector("list", nrow(combs))
  responses <- names(symbol_names)
  M <- lengths(symbol_names)
  C <- length(responses)
  for (i in seq_row(combs)) {
    W_X_B[[i]] <- stats::setNames(vector("list", C), responses)
    for (y in responses) {
      set(d, j = paste0("lag_", y), value = combs[i, y])
    }
    for (y in responses) {
      W_X_B[[i]][[y]] <- model_matrix_emission_formula(
        emission_formula[[y]], d, n_sequences, 1L, 
        n_states, M[y], id_var, time_var, rep(1L, n_sequences), 
        X_mean = attr(X_B[[y]], "X_mean"), check = FALSE, 
        scale = TRUE, R_inv = attr(X_B[[y]], "R_inv")
      )
    }
  }
  W_X_B
}