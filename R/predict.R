#' Predict method for non-homogeneous hidden Markov models
#' 
#' @param object A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param newdata Optional data frame which is used for prediction.
#' @param nsim Non-negative integer defining the number of samples from the
#' normal approximation of the model coefficients.
#' @param probs Vector defining the quantiles of interest.
#' @param return_samples Logical indicating whether to return samples or
#' quantiles. Default is `FALSE`.
#' @param ... Ignored.
#' @export
predict.nhmm <- function(
    object, newdata = NULL, nsim = 0, 
    probs = c(0.025, 0.5, 0.975), return_samples = FALSE, ...) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  if (!is.null(newdata)) {
    time <- object$time_variable
    id <- object$id_variable
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
  } else {
    stopifnot_(
      !is.null(object$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    object <- update(object, newdata = newdata)
  }
  
  beta_i_raw <- stan_to_cpp_initial(
    object$coefficients$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    object$coefficients$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    object$coefficients$beta_o_raw,
    1,
    object$n_channels > 1
  )
  X_initial <- t(object$X_initial)
  X_transition <- aperm(object$X_transition, c(3, 1, 2))
  X_emission <- aperm(object$X_emission, c(3, 1, 2))
  out <- list()
  out$pi <- get_pi(beta_i_raw, X_initial, 0)
  out$A <- get_A(beta_s_raw, X_transition, 0)
  out$B <- if (object$n_channels == 1) {
    get_B(beta_o_raw, X_emission, 0)
  } else {
    get_multichannel_B(
      beta_o_raw, 
      X_emission, 
      object$n_states, 
      object$n_channels, 
      object$n_symbols, 
      0, 0)
  }
  if (nsim > 0) {
    samples <- sample_parameters(object, nsim, probs, return_samples)
    if (return_samples) {
      out$samples <- samples
    } else {
      out$quantiles <- samples
    }
  }
  out
}
#' @export
predict.mnhmm <- function(
    object, newdata = NULL, nsim = 0, 
    probs = c(0.025, 0.5, 0.975), return_samples = FALSE, ...) {
  
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  if (!is.null(newdata)) {
    time <- object$time_variable
    id <- object$id_variable
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
  } else {
    stopifnot_(
      !is.null(object$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    object <- update(object, newdata = newdata)
  }
  
  beta_i_raw <- stan_to_cpp_initial(
    object$coefficients$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    object$coefficients$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    object$coefficients$beta_o_raw,
    1,
    object$n_channels > 1
  )
  X_initial <- t(object$X_initial)
  X_transition <- aperm(object$X_transition, c(3, 1, 2))
  X_emission <- aperm(object$X_emission, c(3, 1, 2))
  X_cluster <- t(object$X_cluster)
  out <- list()
  out$pi <- get_pi(beta_i_raw, X_initial, 0)
  out$A <- get_A(beta_s_raw, X_transition, 0)
  out$B <- if (object$n_channels == 1) {
    get_B(beta_o_raw, X_emission, 0)
  } else {
    get_multichannel_B(
      beta_o_raw, 
      X_emission, 
      object$n_states, 
      object$n_channels, 
      object$n_symbols, 
      0, 0)
  }
  out$omega <- get_omega(
    object$coefficients$theta_raw, X_cluster, 0
  )
  if (nsim > 0) {
    samples <- sample_parameters(object, nsim, probs, return_samples)
    if (return_samples) {
      out$samples <- samples
    } else {
      out$quantiles <- samples
    }
  }
  out
}
