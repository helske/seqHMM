#' Average Marginal Effects for Non-homogenous Hidden Markov Models
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing one or two values for `variable`.
#' @param newdata Optional data frame which is used for marginalization.
#' @param ... Ignored.
#' @rdname ame
#' @export
ame <- function(model, variable, values, ...) {
  UseMethod("ame", model)
}
#' @rdname ame
#' @export
ame.nhmm <- function(
    model, variable, values, newdata = NULL, 
    ...) {
  stopifnot_(
    attr(model, "intercept_only") == FALSE,
    "Model does not contain any covariates."
  )
  stopifnot_(
    checkmate::test_string(x = variable), 
    "Argument {.arg variable} must be a single character string."
  )
  stopifnot_(
    length(values) == 2, 
    "Argument {.arg values} should contain two values for 
    variable {.var variable}.")
  time <- model$time_variable
  id <- model$id_variable
  if (!is.null(newdata)) {
    stopifnot_(
      is.data.frame(newdata), 
      "Argument {.arg newdata} must be a {.cls data.frame} object."
    )
    stopifnot_(
      !is.null(newdata[[id]]), 
      "Can't find grouping variable {.var {id}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[time]]), 
      "Can't find time index variable {.var {time}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[variable]]), 
      "Can't find time variable {.var {variable}} in {.arg newdata}."
    )
  } else {
    stopifnot_(
      !is.null(model$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    newdata <- model$data
  }
  channel <- if (model$n_channels > 1) "channel" else NULL
  
  newdata[[variable]] <- values[1]
  model1 <- update(model, newdata)
  newdata[[variable]] <- values[2]
  model2 <- update(model, newdata)
  
  if (!attr(model, "iv_pi")) {
    X1 <- model1$X_initial[, 1L, drop = FALSE]
    X2 <- model2$X_initial[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_initial
    X2 <- model2$X_initial
  }
  ame_pi <- get_pi_ame(model$gammas$pi, X1, X2, probs)
  if (!attr(model, "iv_A")) {
    X1 <- model1$X_transition[, 1L, drop = FALSE]
    X2 <- model2$X_transition[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_transition
    X2 <- model2$X_transition
  }
  ame_A <- get_A_ame(model$gammas$A, X1, X2, attr(model, "tv_A"), probs)
  if (!attr(model, "iv_B")) {
    X1 <- model1$X_emission[, 1L, drop = FALSE]
    X2 <- model2$X_emission[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_emission
    X2 <- model2$X_emission
  }
  ame_B <- do.call(
    rbind,
    lapply(seq_len(model$n_channels), function(i) {
      get_B_ame(model$gammas$B[[i]], X1, X2, attr(model, "tv_B"), probs)
    })
  )
  out <- list(
    pi = ame_pi,
    A = ame_A,
    B = ame_B
  )
  class(out) <- "amp"
  attr(out, "model") <- "nhmm"
  out
}

#' @rdname ame
#' @export
ame.mnhmm <- function(
    model, variable, values, newdata = NULL, 
    ...) {
  stopifnot_(
    attr(model, "intercept_only") == FALSE,
    "Model does not contain any covariates."
  )
  stopifnot_(
    checkmate::test_string(x = variable), 
    "Argument {.arg variable} must be a single character string."
  )
  stopifnot_(
    length(values) == 2, 
    "Argument {.arg values} should contain two values for 
    variable {.var variable}.")
  time <- model$time_variable
  id <- model$id_variable
  if (!is.null(newdata)) {
    stopifnot_(
      is.data.frame(newdata), 
      "Argument {.arg newdata} must be a {.cls data.frame} object."
    )
    stopifnot_(
      !is.null(newdata[[id]]), 
      "Can't find grouping variable {.var {id}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[time]]), 
      "Can't find time index variable {.var {time}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[variable]]), 
      "Can't find time variable {.var {variable}} in {.arg newdata}."
    )
  } else {
    stopifnot_(
      !is.null(model$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    newdata <- model$data
  }
  
  newdata[[variable]] <- values[1]
  model1 <- update(model, newdata)
  newdata[[variable]] <- values[2]
  model2 <- update(model, newdata)
  D <- model$n_clusters
  
  if (!attr(model, "iv_pi")) {
    X1 <- model1$X_initial[, 1L, drop = FALSE]
    X2 <- model2$X_initial[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_initial
    X2 <- model2$X_initial
  }
  ame_pi <- do.call(
    rbind,
    lapply(seq_len(D), function(i) {
      get_pi_ame(model$gammas$pi[[i]], X1, X2, probs)
    })
  )
  if (!attr(model, "iv_A")) {
    X1 <- model1$X_transition[, 1L, drop = FALSE]
    X2 <- model2$X_transition[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_transition
    X2 <- model2$X_transition
  }
  ame_A <- do.call(
    rbind,
    lapply(seq_len(D), function(i) {
      get_A_ame(model$gammas$A[[i]], X1, X2, attr(model, "tv_A"), probs)
    })
  )
  if (!attr(model, "iv_B")) {
    X1 <- model1$X_emission[, 1L, drop = FALSE]
    X2 <- model2$X_emission[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_emission
    X2 <- model2$X_emission
  }
  ame_B <- do.call(
    rbind,
    lapply(seq_len(D), function(i) {
      do.call(
        rbind,
        lapply(seq_len(model$n_channels), function(j) {
          get_B_ame(model$gammas$B[[i]][[j]], X1, X2, attr(model, "tv_B"), probs)
        })
      )
    })
  )
  out <- list(
    pi = ame_pi,
    A = ame_A,
    B = ame_B
  )
  
  class(out) <- "amp"
  attr(out, "model") <- "mnhmm"
  out
}
