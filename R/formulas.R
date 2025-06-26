#' Get the response variable(s) from a formula(s)
#'
#' @param x A `formula` object or a list of such objects.
#' @param allow_mv A logical value indicating whether to allow multiple 
#' responses in LHS (e.g., c(y, x) ~ 1).
#' @noRd
get_responses <- function(x, allow_mv = TRUE) {
  stopifnot_(
    inherits(x, "formula"), 
    "{.arg emission_formula} must be a {.cls formula} object or a list of 
    {.cls formula} objects."
  )
  stopifnot_(
    identical(length(x), 3L), 
    "{.arg emission_formula} must contain the response variable(s) on the 
    left-hand side of the {.cls formula} object(s)."
  )
  y <- all.vars(x[[2]])
  stopifnot_(
    length(y) == 1L || (length(y) > 1L && allow_mv),
    "{.arg emission_formula} must be a {.cls formula} object with one or more 
    response variables on the left-hand side, or a list of {.cls formula} 
    objects with a single response variable on the LHS of each {.cls formula}"
  )
  y
}

#' Extract lag-terms
#'
#' @param x An RHS of a `formula` object, i.e., a `language` object.
#' @noRd
find_lags <- function(x) {
  if (!is.recursive(x)) {
    return(character(0L))
  }
  if (is.call(x)) {
    if (identical(as.character(x[[1L]]), "lag")) {
      return(deparse1(x))
    } else {
      unlist(lapply(x[-1L], find_lags))
    }
  }
}
#' Extract non-lag-terms
#'
#' @param x  x An RHS of a `formula` object, i.e., a `language` object.
#' @noRd
find_nonlags <- function(x) {
  if (!is.recursive(x)) {
    if (is.name(x)) {
      return(as.character(x))
    }
  }
  if (is.call(x)) {
    if (identical(as.character(x[[1L]]), "lag")) {
      character(0L)
    } else {
      unlist(lapply(x[-1L], find_nonlags))
    }
  }
}
#' Replaces lag(y) terms lag_y in the language object
#' 
#' @param x  x An RHS of a `formula` object, i.e., a `language` object.
#' @noRd
replace_lags <- function(x) {
  if (is.call(x)) {
    if (x[[1]] == as.name("lag") && length(x) == 2 && is.name(x[[2]])) {
      as.name(paste0("lag_", as.character(x[[2]])))
    } else {
      as.call(lapply(x, replace_lags))
    }
  } else {
    x
  }
}
#' Check that a formula for a model is valid and return a modified formula
#' 
#' This function checks that the supplied formula is valid and replaces the 
#' lagged terms of form `lag(variable)` with `lag_variable`.
#' 
#' @param `f` A `formula` object.
#' @param `responses` A character vector of response variables of the model.
#' @param formula_type A character string indicating the type of formula, either
#' `"pi"`, `"A"`, or `"B"`.
#' @param data A `data.table` used to check that variables in formula are present.
#' The data is also updated to include lagged variables.
#' @param id_var Name of the id variable in `data` identifying different 
#' sequences.
#' @noRd
check_formula <- function(f, responses, formula_type, data, id_var) {
  id <- NULL
  if (identical(length(f), 3L)) {
    x <- f[[3L]]
    response <- f[[2L]]
  } else {
    x <- f[[2L]]
    response <- NULL
  }
  lag_terms <- unique(find_lags(x))
  if (formula_type == "pi") {
    stopifnot_(
      length(lag_terms) == 0L || formula_type != "pi",
      "{.arg initial_formula} should not contain any lagged variables."
    )
    stopifnot_(
      is.null(response),
      "{.arg initial_formula} should be a {.cls formula} with empty 
      left-hand side."
    )
    ft <- "initial_formula"
  }
  if (formula_type == "A") {
    stopifnot_(
      is.null(response),
      "{.arg transition_formula} should be a {.cls formula} with empty 
      left-hand side."
    )
    ft <- "transition_formula"
  }
  if (formula_type == "B") {
    ft <- "emission_formula"
  }
  if (formula_type != "pi") {
    vars <- vapply(lag_terms, \(i) sub("lag\\((.*?)\\)", "\\1", i), "")
    z <- vars[which(!(vars %in% responses))]
    stopifnot_(
      length(z) == 0L,
      c(
        "{cli::qty(responses)} Only lagged response variable{?s} {responses} 
        are allowed in {.arg {ft}}. ",
        `x` = "For lagged predictor{?s} {z}, create lagged variable{?s} as 
        {?a/} new column{?s} in the {.arg data}."
      )
    )
  } else {
    vars <- NULL
  }
  
  nonlag_terms <- unique(find_nonlags(x))
  z <- nonlag_terms[which(nonlag_terms %in% responses)]
  stopifnot_(
    length(z) == 0L,
    c("RHS of model formulas should not contain non-lagged response variables.",
      `x` = "{cli::qty(z)} Found {?a/} response variable{?s} {z} in the `{ft}`. "
    )
  )
  stopifnot_(
    all(nonlag_terms %in% colnames(data)), 
    "Not all variables defined in `{ft}` are present in {.arg data} ."
  )
  f <- stats::as.formula(replace_lags(f), env = environment(f))
  if (length(vars) > 0L) {
    for (y in vars) {
      lag_obs <- paste0("lag_", y)
      if (is.null(data[[lag_obs]])) {
        y1 <- data[[y]][1] # this value is not used anywhere
        data[, lag_obs := shift(y, type = "lag", fill = y1), by = id, 
             env = list(id = id_var, y = y, lag_obs = lag_obs, y1 = y1),
             showProgress = FALSE]
      }
    }
    attr(f, "responses") <- vars
  }
  f
}
