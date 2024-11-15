#' Average Marginal Effects for NHMM Responses
#' 
#' The function `ame_obs` computes the average marginal effect (AME) of the 
#' model covariate \eqn{X} at time t on the current and future responses by 
#' marginalizing over the sequences and latent states. Under the assumption of 
#' no unobserved confounding (i.e., there are no unobserved variables that 
#' influence the covariate \eqn{X} and the outcome \eqn{Y}), these can be 
#' regarded as the causal effects. In case `values` argument is a single value 
#' \eqn{x}, the function returns the interventional distributions
#' \deqn{P(y_{t+k} | do(X_t = x))}
#' and in a case `values` contains two values \eqn{x} and \eqn{w} a shift in 
#' interventional distributions, i.e.,
#' \deqn{P(y_{t+k} | do(X_t = x)) - P(y_{t+k} | do(X_t = w))}.
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing one or two values for `variable`. 
#' See details.
#' @param start_time Time(s) of intervention. Either a scalar or vector. 
#' Intervention is applied to all provided time points. 
#' @param newdata Optional data frame which is used for marginalization.
#' @param probs Quantiles of interest of average marginal effect.
#' @param ... Ignored.
#' @rdname ame_obs
#' @export
ame_obs <- function(model, variable, values, start_time, ...) {
  UseMethod("ame_obs", model)
}
#' @rdname ame_obs
#' @export
ame_obs.nhmm <- function(
    model, variable, values, start_time, newdata = NULL, probs = c(0.05, 0.95),
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
  stopifnot_(
    !is.null(model$boot),
    paste0(
      "Model does not contain bootstrap samples of coefficients. ",
      "Run {.fn bootstrap_coefs} first."
    )
  )
  newdata[[variable]] <- values[1]
  model1 <- update(model, newdata)
  newdata[[variable]] <- values[2]
  model2 <- update(model, newdata)
  C <- model$n_channels
  if (C == 1L) {
    times <- colnames(model$observations)
    symbol_names <- list(model$symbol_names)
  } else {
    times <- colnames(model$observations[[1]])
    symbol_names <- model$symbol_names
  }
  stop("WIP")
  if (model$n_channels == 1) {
    
    obs <- create_obsArray(model)[1L, , ]
    out1 <- state_obs_probs_nhmm_singlechannel( 
      model1$etas$pi, model1$X_pi, model1$etas$A, model1$X_A, 
      model1$etas$B, model1$X_B, obs, model1$sequence_lengths, 
      attr(model1$X_pi, "icpt_only"), attr(model1$X_A, "icpt_only"), 
      attr(model1$X_B, "icpt_only"), attr(model1$X_A, "iv"), 
      attr(model1$X_B, "iv"), attr(model1$X_A, "tv"), attr(model1$X_B, "tv"),
      start = start_time)
    out2 <- state_obs_probs_nhmm_singlechannel( 
      model2$etas$pi, model2$X_pi, model2$etas$A, model2$X_A, 
      model2$etas$B, model2$X_B, obs, model2$sequence_lengths, 
      attr(model2$X_pi, "icpt_only"), attr(model2$X_A, "icpt_only"), 
      attr(model2$X_B, "icpt_only"), attr(model2$X_A, "iv"), 
      attr(model2$X_B, "iv"), attr(model2$X_A, "tv"), attr(model2$X_B, "tv"),
      start = start_time)
  }
  
  class(out) <- "ame_obs"
  attr(out, "model") <- "nhmm"
  out
}

#' @rdname ame_obs
#' @export
ame_obs.mnhmm <- function(
    model, variable, values, start_time, newdata = NULL, probs = c(0.05, 0.95),
    ...) {
  
  stop("Not yet implemented")
  class(out) <- "ame_obs"
  attr(out, "model") <- "mnhmm"
  out
}
