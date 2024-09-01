#' Average Marginal Effects of Covariates of Non-homogenous Hidden Markov Models
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing two values for `variable`.
#' @param newdata Optional data frame which is used for marginalization.
#' @param nsim Non-negative integer defining the number of samples from the 
#' normal approximation of the model parameters used in 
#' computing the approximate quantiles of the estimates. If `0`, only point 
#' estimates are returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.975)`.
#' @export
average_marginal_effect <- function(
    model, variable, values, newdata = NULL, 
    nsim = 0, 
    probs = c(0.025, 0.975)) {
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  stopifnot_(
    checkmate::test_string(x = variable), 
    "Argument {.arg variable} must be a single character string."
  )
  stopifnot_(
    length(values) != 2, 
    "Argument {.arg values} should contain two values for 
    variable {.var variable}.")
  if (is.null(newdata)) {
    time <- model$time_variable
    id <- model$id_variable
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
    stopifnot(
      !is.null(object$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    newdata <- model$data
  }
  
  beta_i_raw <- stan_to_cpp_initial(
    model$estimation_results$parameters$beta_i_raw
  )
  beta_s_raw <- stan_to_cpp_transition(
    model$estimation_results$parameters$beta_s_raw
  )
  beta_o_raw <- stan_to_cpp_emission(
    model$estimation_results$parameters$beta_o_raw,
    1,
    C > 1
  )
  newdata[[variable]] <- values[1]
  model <- update(model, newdata = newdata)
  X_initial1 <- t(model$X_initial)
  X_transition1 <- aperm(model$X_transition, c(3, 1, 2))
  X_emission1 <- aperm(model$X_emission, c(3, 1, 2))
  newdata[[variable]] <- values[2]
  model <- update(model, newdata = newdata)
  X_initial2 <- t(model$X_initial)
  X_transition2 <- aperm(model$X_transition, c(3, 1, 2))
  X_emission2 <- aperm(model$X_emission, c(3, 1, 2))
  
  ame_pi <- get_pi(beta_i_raw, X_initial1, 0) - 
    get_pi(beta_i_raw, X_initial2, 0)
  ame_A <- get_A(beta_s_raw, X_transition1, 0) - 
    get_A(beta_s_raw, X_transition2, 0)
  ame_B <- if (model$n_channels == 1) {
    get_B(beta_o_raw, X_emission1, 0) - get_B(beta_o_raw, X_emission2, 0)
  } else {
    get_multichannel_B(beta_o_raw, X_emission1, S, C, M, 0) -
      get_multichannel_B(beta_o_raw, X_emission2, S, C, M, 0) 
  }
  browser()
  if (nsim > 0) {
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    chol_precision <- chol(-model$estimation$hessian)
    U <- backsolve(chol_precision, diag(ncol(chol_precision)))
    x <- matrix(rnorm(nsim * ncol(U)), nrow = nsim) %*% U
    x <- t(sweep(x, 2, c(beta_i_raw, beta_s_raw, beta_o_raw), "+"))
    p_i <- length(beta_i_raw)
    p_s <- length(beta_s_raw)
    p_o <- length(beta_o_raw)
    
    pi_samples <- apply(
      x[seq_len(p_i), ], 2, function(z) {
        z <- array(z, dim = dim(beta_i_raw))
        z <- stan_to_cpp_initial(z)
        get_pi(z, X_initial1) - get_pi(z, X_initial2)
      }
    )
    A_samples <- apply(
      x[p_i + seq_len(p_s), ], 2, function(z) {
        z <- array(z, dim = dim(beta_s_raw))
        z <- stan_to_cpp_transition(z)
        unlist(get_A(z, X_transition1)) - 
          unlist(get_A(z, X_transition2))
      }
    )
    B_samples <- apply(
      x[p_i + p_s + seq_len(p_o), ], 2, function(z) {
        z <- array(z, dim = dim(beta_o_raw))
        z <- stan_to_cpp_emission(z, 1, model$n_channels > 1)
        unlist(get_B(aperm(z, c(2, 3, 1)), X_emission1)) -
          unlist(get_B(aperm(z, c(2, 3, 1)), X_emission2))
      }
    )
    
    quantiles <- fast_quantiles(samples, probs)
    for(i in seq_along(probs)) {
      transition_probs[paste0("q", 100 * probs[i])] <- quantiles[, i]
    }
  }
  class(out) <- "ame"
}
