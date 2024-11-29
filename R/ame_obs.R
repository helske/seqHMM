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
    model, variable, values, start_time, newdata = NULL, probs, ...) {
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
  if (!missing(probs)) {
    return_quantiles <- TRUE
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first in order to compute quantiles."
      )
    )
  } else {
    return_quantiles <- FALSE
    # dummy stuff for C++
    model$boot <- list(
      gamma_pi = list(model$gammas$pi),
      gamma_A = list(model$gammas$A),
      gamma_B = list(model$gammas$B)
    )
    probs <- 0.5
  }
  
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
  newdata[[variable]][newdata[[time]] >= start_time] <- values[1]
  X1 <- update(model, newdata)[c("X_pi", "X_A", "X_B")]
  newdata[[variable]][newdata[[time]] >= start_time] <- values[2]
  X2 <- update(model, newdata)[c("X_pi", "X_A", "X_B")]
  C <- model$n_channels
  if (C == 1L) {
    times <- as.numeric(colnames(model$observations))
    symbol_names <- list(model$symbol_names)
    obs <- create_obsArray(model)[1L, , ]
    out <- ame_obs_nhmm_singlechannel( 
      model$etas$pi, model$etas$A, model$etas$B, obs, 
      model$sequence_lengths, 
      attr(X1$X_pi, "icpt_only"), attr(X1$X_A, "icpt_only"), 
      attr(X1$X_B, "icpt_only"), attr(X1$X_A, "iv"), 
      attr(X1$X_B, "iv"), attr(X1$X_A, "tv"), attr(X1$X_B, "tv"),
      X1$X_pi, X1$X_A, X1$X_B, 
      attr(X2$X_pi, "icpt_only"), attr(X2$X_A, "icpt_only"), 
      attr(X2$X_B, "icpt_only"), attr(X2$X_A, "iv"), 
      attr(X2$X_B, "iv"), attr(X2$X_A, "tv"), attr(X2$X_B, "tv"), 
      X2$X_pi, X2$X_A, X2$X_B,
      model$boot$gamma_pi, model$boot$gamma_A, model$boot$gamma_B, start_time, 
      probs
    )
    d <- data.frame(
      observation = model$symbol_names,
      time = rep(as.numeric(colnames(model$observations)), each = model$n_symbols),
      estimate = c(out$point_estimate)
    )
    if (return_quantiles) {
      for(i in seq_along(probs)) {
        d[paste0("q", 100 * probs[i])] <- c(out$quantiles[, , i])
      }
    }
  } else {
    times <- as.numeric(colnames(model$observations[[1]]))
    symbol_names <- model$symbol_names
    obs <- create_obsArray(model)
    out <- ame_obs_nhmm_multichannel( 
      model$etas$pi, model$etas$A, model$etas$B, obs, 
      model$sequence_lengths, 
      attr(X1$X_pi, "icpt_only"), attr(X1$X_A, "icpt_only"), 
      attr(X1$X_B, "icpt_only"), attr(X1$X_A, "iv"), 
      attr(X1$X_B, "iv"), attr(X1$X_A, "tv"), attr(X1$X_B, "tv"),
      X1$X_pi, X1$X_A, X1$X_B, 
      attr(X2$X_pi, "icpt_only"), attr(X2$X_A, "icpt_only"), 
      attr(X2$X_B, "icpt_only"), attr(X2$X_A, "iv"), 
      attr(X2$X_B, "iv"), attr(X2$X_A, "tv"), attr(X2$X_B, "tv"), 
      X2$X_pi, X2$X_A, X2$X_B,
      model$boot$gamma_pi, model$boot$gamma_A, model$boot$gamma_B, start_time, 
      probs
    )
    d <- data.frame(
      observation = model$symbol_names,
      time = rep(as.numeric(colnames(model$observations)), each = model$n_symbols),
      estimate = c(out$point_estimate)
    )
    if (return_quantiles) {
      for(i in seq_along(probs)) {
        d[paste0("q", 100 * probs[i])] <- c(out$quantiles[, , i])
      }
    }
  }
  d[d[[time]] >= start_time, ]
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
