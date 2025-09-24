#' Functions for computing the marginal probabilities for `get_marginals`
#' @noRd
compute_z_marginals <- function(model, id_time, pp, cond, cluster = NULL) {
  # avoid CRAN check warnings due to NSE
  probability <- cols <- NULL
  x <- model$data[, cols, env = list(cols = as.list(c(id_time, cond)))]
  x <- pp[x, on = id_time, nomatch = 0L]
  cond <- c(cond, cluster, "state")
  x[, list(probability = mean(probability)), by = cond,
    showProgress = FALSE]
}
compute_A_marginals <- function(model, id_time, pp, cond, cluster = NULL) {
  # avoid CRAN check warnings due to NSE
  probability <- i.probability <- state_prob <- cols <- NULL
  A <- get_transition_probs(model)
  x <- model$data[, cols, env = list(cols = as.list(c(id_time, cond)))]
  x <- A[x, on = id_time, nomatch = 0L]
  x[pp, state_prob := i.probability, on = c(id_time, state_from = "state")]
  cond <- c(cond, cluster, "state_from", "state_to")
  x[, list(probability = sum(probability * state_prob) / sum(state_prob)), 
    by = cond, showProgress = FALSE]
}
compute_y_and_B_marginals <- function(model, id_time, pp, cond, cluster = NULL) {
  # avoid CRAN check warnings due to NSE
  probability <- i.probability <- state_prob <- cols <- NULL
  B <- get_emission_probs(model)
  x <- model$data[, cols, env = list(cols = as.list(c(id_time, cond)))]
  p_B <- p_y <- vector("list", model$n_channels)
  for (i in seq_len(model$n_channels)) {
    x_i <- B[[i]][x, on = id_time, nomatch = 0L]
    x_i[pp, state_prob := i.probability, on = c(id_time, "state")]
    cond_i <- c(cond, model$responses[i])
    p_y[[i]] <- x_i[, list(probability = sum(probability * state_prob) / 
                             sum(state_prob)), by = cond_i,
                    showProgress = FALSE]
    cond_i <- c(cluster, "state", cond_i)
    p_B[[i]] <- x_i[, list(probability = sum(probability * state_prob) / 
                             sum(state_prob)), by = cond_i,
                    showProgress = FALSE]
  }
  list(B = p_B, y = p_y)
}

#' Compute the Marginal Probabilities from NHMMs
#' 
#' `get_marginals` returns the marginal state, response, transition, and emission
#' probabilities, optionally per grouping defined by `condition`. By default, 
#' the marginalization weights sequences by the corresponding posterior 
#' probabilities of the latent states, i.e., conditional probabilities of the 
#' latent states given all data (`weighting = "posterior"`). If 
#' `weighting = "forward"`, marginalization is based on forward probabilities, 
#' i.e. state probabilities given data up to that point which allows you to 
#' compute, for example, state marginals of form 
#' \eqn{P(state_t | data_1, \ldots, data_t)} (whereas in posterior probability 
#' weighting the conditioning is on \eqn{data_1,\ldots,data_T}. 
#' If `weighting = "none"`, all individuals and time points are treated equally, 
#' without accounting for the probability that individual is at particular 
#' state at particular time.
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `NULL`, in which case no quantiles are computed. The quantiles are based on 
#' bootstrap samples of coefficients, stored in `object$boot`.
#' @param condition An optional vector of variable names used for conditional 
#' marginal probabilities. Default is `NULL`, in which case marginalization is 
#' done over all variables, so that for example marginal emission probabilities 
#' are computed over all individuals and time points.
#' @param newdata An optional data frame containing the new data to be used in 
#' computing the probabilities.
#' @param type A character vector defining the marginal probabilities of 
#' interest. Can be one or multiple of `"state"`, `"response"`, `"transition"`, 
#' and `"emission"`. Default is to compute all of these.
#' @param weighting A character string defining the type of weighting used in 
#' marginalization. One of `"posterior"` , `"forward"`, `"none"`. See details.
#' @export
get_marginals <- function(model, probs = NULL, condition = NULL, 
                          newdata = NULL, 
                          type = c("state", "response", "transition", "emission"),
                          weighting = c("posterior", "forward", "none")) {
  
  # avoid CRAN check warnings due to NSE
  time <- probability <- NULL
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
  out_state <- out_A <- out_obs <- out_B <- NULL
  if (inherits(model, "mnhmm")) {
    cluster <- "cluster"
  } else {
    cluster <- NULL
  }
  if (compute_z <- "state" %in% type) {
    out_state <- compute_z_marginals(model, id_time, pp, condition, cluster)
  }
  if (compute_A <- "transition" %in% type) {
    out_A <- compute_A_marginals(model, id_time, pp, condition, cluster)
  }
  compute_y <- compute_B <- FALSE
  if ("response" %in% type || "emission" %in% type) {
    out_y_B <- compute_y_and_B_marginals(model, id_time, pp, condition, cluster)
    if (compute_y <- "response" %in% type) {
      out_obs <- stats::setNames(out_y_B$y, model$responses)
    }
    if (compute_B <- "emission" %in% type) {
      out_B <- stats::setNames(out_y_B$B, model$responses)
    }
  }
  if (!is.null(probs)) {
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    nsim <- length(model$boot$gamma_pi)
    boot_state <- matrix(0, nrow(out_state), nsim * compute_z)
    boot_A <- matrix(0, nrow(out_A), nsim * compute_A)
    boot_obs <- vector("list", model$n_channels)
    boot_B <- vector("list", model$n_channels)
    for (i in seq_len(model$n_channels)) {
      boot_B[[i]] <- matrix(0, nrow(out_B[[i]]), nsim * compute_B)
      boot_obs[[i]] <- matrix(0, nrow(out_obs[[i]]), nsim * compute_y)
    }
    tQs <- t(create_Q(model$n_states))
    tQm <- lapply(model$n_symbols, \(i) t(create_Q(i)))
    for (i in seq_len(nsim)) {
      model$gammas$gamma_pi[] <- model$boot$gamma_pi[[i]]
      model$gammas$gamma_A[] <- model$boot$gamma_A[[i]]
      model$gammas$gamma_B <- model$boot$gamma_B[[i]]
      if (inherits(model, "mnhmm")) {
        model$gammas$gamma_omega[] <- model$boot$gamma_omega[[i]]
      }
      if (weighting == "posterior") {
        pp <- posterior_probs(model)
      }
      if (weighting == "forward") {
        pp <- forward_backward(model, forward_only = TRUE)
        setnames(pp, "log_alpha", "probability")
      }
      if (compute_z) {
        boot_state[, i] <- compute_z_marginals(
          model, id_time, pp, condition, cluster
        )$probability
      }
      if (compute_A) {
        boot_A[, i] <- compute_A_marginals(
          model, id_time, pp, condition, cluster
        )$probability
      }
      if (compute_y || compute_B) {
        out_y_B <- compute_y_and_B_marginals(
          model, id_time, pp, condition, cluster
        )
        if (compute_y) {
          for (j in seq_len(model$n_channels)) {
            boot_obs[[j]][, i] <- out_y_B$y[[j]]$probability
          }
        }
        if (compute_B) {
          for (j in seq_len(model$n_channels)) {
            boot_B[[j]][, i] <- out_y_B$B[[j]]$probability
          }
        }
      }
    }
    if (compute_z) {
      qs <- t(apply(boot_state, 1, quantileq, probs = probs))
      out_state <- cbind(out_state, qs)
    }
    if (compute_A) {
      qs <- t(apply(boot_A, 1, quantileq, probs = probs))
      out_A <- cbind(out_A, qs)
    }
    if (compute_y) {
      for (i in seq_len(model$n_channels)) {
        qs <- t(apply(boot_obs[[i]], 1, quantileq, probs = probs))
        out_obs[[i]] <- cbind(out_obs[[i]], qs)
      }
    }
    if (compute_B) {
      for (i in seq_len(model$n_channels)) {
        qs <- t(apply(boot_B[[i]], 1, quantileq, probs = probs))
        out_B[[i]] <- cbind(out_B[[i]], qs)
      }
    }
  }
  list(states = out_state, responses = out_obs, transitions = out_A, emissions = out_B)
}
