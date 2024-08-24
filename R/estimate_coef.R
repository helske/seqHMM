#' Estimate Only the Regression Coefficients of Mixture Hidden Markov Models
#'
#' Function `estimate_coef` estimates the regression coefficients of 
#' mixture hidden Markov models of class `mhmm` and its restricted variants 
#' while keeping other parameters fixed.
#'
#' @export
#' @param model An object of class `mhmm`.
#' @param threads Number of threads to use in parallel computing. 
#' The default is 1.
estimate_coef <- function(model, threads = 1) {
  stopifnot(
    inherits(model, "mhmm"),
    "{.arg model} must be a {.cls mhmm} object."
  )
  check_positive_integer(threads)
  df <- attr(model, "df")
  nobs <- attr(model, "nobs")
  original_model <- model
  model <- .combine_models(model)
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)

  em.con <- list(print_level = 0, maxeval = 1000, reltol = 1e-10)

  res <- estimate_coefs(
    model$transition_probs, emissionArray, model$initial_probs, obsArray,
    model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters,
    em.con$maxeval, em.con$reltol, em.con$print_level, threads
  )
  if (res$error != 0) {
    err_msg <- switch(res$error,
      "Scaling factors contain non-finite values.",
      "Backward probabilities contain non-finite values.",
      "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
      "Estimation of coefficients of covariates failed due to singular Hessian.",
      "Estimation of coefficients of covariates failed due to non-finite cluster probabilities.",
      "Non-finite log-likelihood"
    )
    stop_(paste("EM algorithm failed:", err_msg))
  }
  original_model$coefficients[] <- res$coefficients
  original_model
}
