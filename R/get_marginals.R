compute_z_marginals <- function(model, id_time, pp, cond) {
  probability <- cols <- NULL
  x <- model$data[, cols, env = list(cols = as.list(c(id_time, cond)))]
  x <- pp[x, on = id_time, nomatch = 0L]
  cond <- c(cond, "state")
  x[, list(probability = mean(probability)), by = cond]
}

compute_y_and_B_marginals <- function(model, id_time, pp, cond) {
  probability <- i.probability <- state_prob <- cols <- NULL
  B <- get_emission_probs(model)
  x <- model$data[, cols, env = list(cols = as.list(c(id_time, cond)))]
  p_B <- p_y <- vector("list", model$n_channels)
  for (i in seq_len(model$n_channels)) {
    x_i <- B[[i]][x, on = id_time, nomatch = 0L]
    x_i[pp, state_prob := i.probability, on = c(id_time, "state")]
    cond_i <- c(cond, model$responses[i])
    p_y[[i]] <- x_i[, list(probability = sum(probability * state_prob) / sum(state_prob)), by = cond_i]
    cond_i <- c("state", cond_i)
    p_B[[i]] <- x_i[, list(probability = sum(probability * state_prob) / sum(state_prob)), by = cond_i]
  }
  list(B = p_B, y = p_y)
}

compute_A_marginals <- function(model, id_time, pp, cond) {
  probability <- i.probability <- state_prob <- cols <- NULL
  A <- get_transition_probs(model)
  x <- model$data[, cols, env = list(cols = as.list(c(id_time, cond)))]
  x <- A[x, on = id_time, nomatch = 0L]
  x[pp, state_prob := i.probability, on = c(id_time, state_from = "state")]
  cond <- c(cond, "state_from", "state_to")
  x[, list(probability = sum(probability * state_prob) / sum(state_prob)), by = cond]
}

#' Compute the Marginal Probabilities from NHMMs
#' 
#' `get_marginals` returns the marginal state, response, transition, and emission
#' probabilities, optionally per grouping defined by `condition`. By default, 
#' the marginalization weights sequences by the corresponding posterior 
#' probabilities of the latent states (`weighting = "newdata"`). If 
#' `weighting = "original"` or `newdata` is `NULL`, the posterior probabilities 
#' are computed using the data in `model$data`. If `weighting = "none"`, all 
#' individuals and time points are treated equally, without accounting for the 
#' probability that individual is at particular state at particular time. The 
#' option `weighting = "original"` is mostly useful if you are interested in 
#' obtaining one-step marginals with new data, e.g., 
#' \eqn{P(state_t | newdata_t, state_{t-1}, originaldata_{t-1},\ldots, original_data_1)}.
#' 
#' @param model An object of class `nhmm`.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `NULL`, in which case no quantiles are computed. The quantiles are based on 
#' bootstrap samples of coefficients, stored in `object$boot`.
#' @param condition An optional vector of variable names used for conditional 
#' marginal probabilities. Default is `NULL`, in which case no marginalization is 
#' done.
#' @param newdata An optional data frame containing the new data to be used in 
#' computing the probabilities.
#' @param weighting A character string defining the type of weighting used in 
#' marginalization. See details.
#' @export
get_marginals <- function(model, probs = NULL, condition = NULL, 
                          newdata = NULL, 
                          weighting = c("posterior", "forward", "none")) {
  
  # avoid CRAN check warnings due to NSE
  time <- NULL
  stopifnot_(
    inherits(model, c("nhmm", "mnhmm")),
    "Argument {.arg model} must be an object of class {.cls nhmm} or {.cls mnhmm}."
  )
  weighting <- try(
    match.arg(weighting, c("posterior", "forward", "none")), 
    silent = TRUE
  )
  stopifnot_(
    !inherits(weighting, "try-error"),
    "Argument {.arg weighting} must be {.val posterior}, {.val forward}, or 
    {.val none}."
  )
  if (!is.null(condition)) {
    if(is.null(newdata)) {
      stopifnot_(
        all(condition %in% colnames(model$data)), 
        "Not all variables defined in {.arg condition} are present in {.arg model$data} ."
      )
    } else {
      stopifnot_(
        all(condition %in% colnames(newdata)), 
        "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
      )
    }
  }
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  if (weighting == "posterior") {
    pp <- posterior_probs(model)
  }
  if (weighting == "forward") {
    pp <- forward_backward(model, forward_only = TRUE)
    setnames(pp, "log_alpha", "probability")
  }
  if (weighting == "none") {
    pp[, probability := 1L]
  }
  
  id_time <- c(model$id_variable, model$time_variable)
  out_state <- compute_z_marginals(model, id_time, pp, condition)
  out_A <- compute_A_marginals(model, id_time, pp, condition)
  out_y_B <- compute_y_and_B_marginals(model, id_time, pp, condition)
  out_obs <- stats::setNames(out_y_B$y, model$responses)
  out_B <- stats::setNames(out_y_B$B, model$responses)
  
  if (!is.null(probs)) {
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    nsim <- length(model$boot$gamma_pi)
    boot_state <- matrix(0, nrow(out_state), nsim)
    boot_A <- matrix(0, nrow(out_A), nsim)
    boot_obs <- matrix(0, nrow(out_obs), nsim)
    boot_B <- vector("list", model$n_channels)
    for (i in seq_len(model$n_channels)) {
      boot_B[[i]] <- matrix(0, nrow(out_B[[i]]), nsim)
    }
    tQs <- t(create_Q(model$n_states))
    tQm <- t(create_Q(model$n_states))
    for (i in seq_len(nsim)) {
      model$gammas$pi[] <- model$boot$gamma_pi[[i]]
      model$gammas$A[] <- model$boot$gamma_A[[i]]
      model$gammas$B[] <- model$boot$gamma_B[[i]] #Fix for multichannel TODO
      model$etas$pi[] <- tQs %*% model$gammas$pi[]
      for (s in seq_len(model$n_states)) {
        model$etas$A[, , s] <- tQs %*% model$gammas$A[, , s]
        model$etas$B[, , s] <- tQm %*% model$gammas$B[, , s]
      }
      if (weighting == "posterior") {
        pp <- posterior_probs(model)
      }
      if (weighting == "forward") {
        pp <- forward_backward(model, forward_only = TRUE)
        setnames(pp, "log_alpha", "probability")
      }
      boot_state[, i] <- compute_z_marginals(
        model, id_time, pp, condition)$probability
      boot_A[, i] <- compute_A_marginals(model, id_time, pp, condition)$probability
      out_y_B <- compute_y_and_B_marginals(model, id_time, pp, condition)
      boot_obs[, i] <- out_y_B$y$probability
      boot_B[, i] <- out_y_B$B$probability
    }
    qs <- t(apply(boot_state, 1, quantileq, probs = probs))
    out_state <- cbind(out_state, qs)
    qs <- t(apply(boot_A, 1, quantileq, probs = probs))
    out_A <- cbind(out_A, qs)
    qs <- t(apply(boot_obs, 1, quantileq, probs = probs))
    out_obs <- cbind(out_obs, qs)
    qs <- t(apply(boot_B, 1, quantileq, probs = probs))
    out_B <- cbind(out_B, qs)
  }
  list(states = out_state, responses = out_obs, transitions = out_A, emissions = out_B)
}