#' Compute marginal and conditional probabilities from joint distributions 
#' obtained from predict based on forward predictions
#' @noRd
predict_y_B <- function(x, newdata, response, cond) {
  # avoid CRAN check warnings due to NSE
  estimate <- probability <- NULL
  set(newdata, j = "estimate", value = x)
  cond_all <- c(cond, "state", response)
  # Compute joint distribution P(Y_t, Z_t)
  d <- newdata[, list(probability = mean(estimate)), by = cond_all, 
               showProgress = FALSE]
  setorderv(d, cond_all)
  out_y <- d[, list(probability = sum(probability)), 
             by = c(response, cond), 
             showProgress = FALSE]
  out_B <- d[, "probability" := list(probability / sum(probability)), 
             by = c("state", cond),
             showProgress = FALSE]
  
  list(y = out_y, B = out_B)
}
predict_z_A <- function(x, newdata, cond) {
  # avoid CRAN check warnings due to NSE
  estimate <- probability <- NULL
  set(newdata, j = "estimate", value = x)
  cond_all <- c(cond, "state_from", "state_to")
  # Compute joint distribution P(Z_t-1, Z_t)
  d <- newdata[, list(probability = mean(estimate)), by = cond_all,
               showProgress = FALSE]
  setorderv(d, cond_all)
  out_z <- d[, list(probability = sum(probability)), by = c("state_to", cond),
             showProgress = FALSE]
  setnames(out_z, "state_to", "state")
  out_A <- d[, "probability" := list(probability / sum(probability)), 
             by = c("state_from", cond),
             showProgress = FALSE]
  list(z = out_z, A = out_A)
}

.predict_nhmm <- function(
    dy, dz, condition, obs, object, type_y, type_B, type_z, type_A) {
  
  if (inherits(object, "fanhmm")) {
    out <- Rcpp_predict_fanhmm(
      obs, object$sequence_lengths, object$n_symbols, 
      object$X_pi, object$X_A, object$X_B, 
      io(object$X_pi), io(object$X_A), io(object$X_B), 
      iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
      object$gammas$gamma_pi, object$gammas$gamma_A, object$gammas$gamma_B,
      object$prior_obs, object$W_X_B, object$W_A, object$W_B
    )
  } else {
    out <- Rcpp_predict_nhmm(
      obs, object$sequence_lengths, object$n_symbols, 
      object$X_pi, object$X_A, object$X_B, 
      io(object$X_pi), io(object$X_A), io(object$X_B), 
      iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
      object$gammas$gamma_pi, object$gammas$gamma_A, object$gammas$gamma_B
    )
  }
  results_y <- results_B <- stats::setNames(
    vector("list", object$n_channels), object$responses
  )
  results_z <- results_A <- NULL
  if (type_y || type_B) {
    for (y in seq_along(object$responses)) {
      d <- predict_y_B(
        unlist(out$observations[y, ,]), dy[[y]], object$responses[y], 
        condition
      )
      if (type_y) {
        results_y[[y]] <- d$y  
      }
      if (type_B) {
        results_B[[y]] <- d$B
      }
    }
  }
  if (type_z || type_A) {
    d <- predict_z_A(unlist(out$states), dz, condition)
    if (type_z) {
      results_z <- d$z
    }
    if (type_A) {
      results_A <- d$A
    }
  }
  list(y = results_y, B = results_B, z = results_z, A = results_A)
}
.predict_mnhmm <- function(
    dy, dz, condition, obs, object, type_y, type_B, type_z, type_A) {
  # avoid CRAN check warnings due to NSE
  probability <- NULL
  if (inherits(object, "fanhmm")) {
    out <- Rcpp_predict_mfanhmm(
      obs, object$sequence_lengths, object$n_symbols, 
      object$X_pi, object$X_A, object$X_B, object$X_omega,
      io(object$X_pi), io(object$X_A), io(object$X_B), io(object$X_omega),
      iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
      object$gammas$gamma_pi, object$gammas$gamma_A, object$gammas$gamma_B,
      object$gammas$gamma_omega,
      object$prior_obs, object$W_X_B, object$W_A, object$W_B
    )
  } else {
    out <- Rcpp_predict_mnhmm(
      obs, object$sequence_lengths, object$n_symbols, 
      object$X_pi, object$X_A, object$X_B, object$X_omega,
      io(object$X_pi), io(object$X_A), io(object$X_B), io(object$X_omega),
      iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
      object$gammas$gamma_pi, object$gammas$gamma_A, object$gammas$gamma_B,
      object$gammas$gamma_omega
    )
  }
  results_z <- results_A <- results_y <- results_B <- stats::setNames(
    vector("list", object$n_clusters), object$cluster_names
  )
  if (type_y || type_B) {
    for (i in seq_along(object$cluster_names)) {
      results_y[[i]] <- results_B[[i]] <- stats::setNames(
        vector("list", object$n_channels), object$responses
      )
      states <- object$state_names[[i]]
      for (y in seq_along(object$responses)) {
        M <- object$n_symbols[y]
        d <- predict_y_B(
          unlist(out$observations[i, y,]), dy[[y]], object$responses[y], 
          condition
        )
        if (type_y) {
          results_y[[i]][[y]] <- d$y[, probability := probability / sum(probability)]
        }
        if (type_B) {
          results_B[[i]][[y]] <- d$B
          set(results_B[[i]][[y]], j = "state", value = rep(states, each = M))
        }
      }
    }
  }
  if (type_z || type_A) {
    for (i in seq_along(object$cluster_names)) {
      d <- predict_z_A(unlist(out$states[i, ,]), dz, condition)
      if (type_z) {
        results_z[[i]] <- d$z[, probability := probability / sum(probability)] 
      }
      if (type_A) {
        S <- length(states)
        results_A[[i]] <- d$A
        set(results_A[[i]], j = "state_from", value = rep(states, each = S))
        set(results_A[[i]], j = "state_to", value = rep(states, S))
      }
    }
  }
  list(y = results_y, B = results_B, z = results_z, A = results_A)
}
prob_diff <- function(d1, d2, responses = NULL, clusters = NULL) {
  if (is.null(clusters)) {
    if (is.null(responses)) {
      set(d1, j = "probability", 
          value = d1[["probability"]] - d2[["probability"]])
    } else {
      for (y in responses) {
        set(d1[[y]], j = "probability", 
            value = d1[[y]][["probability"]] - d2[[y]][["probability"]])
      }
    }
  } else {
    for (i in seq_along(clusters)) {
      if (is.null(responses)) {
        set(d1[[i]], j = "probability", 
            value = d1[[i]][["probability"]] - d2[[i]][["probability"]])
      } else {
        for (y in responses) {
          set(d1[[y]], j = "probability", 
              value = d1[[i]][[y]][["probability"]] - d2[[i]][[y]][["probability"]])
        }
      }
    }
  }
}


#' Predictions from Non-homogeneous Hidden Markov Models
#'
#' This function computes the marginal forward predictions for NHMMs and 
#' MNHMMs, where the marginalization is (by default) over individuals and time 
#' points, weighted by the latent state probabilities.

#' @param object An object of class `nhmm` or `mnhmm`.
#' @param newdata A data frame used for computing the predictions.
#' @param newdata2 An optional data frame for predictions, in which case the 
#' estimates are differences between predictions using `newdata` and `newdata2`.
#' @param condition An optional vector of variable names used for conditional 
#' predictions.
#' @param type A character vector defining the marginal predictions of 
#' interest. Can be one or multiple of `"state"`, `"response"`, `"transition"`, 
#' and `"emission"`. Default is to compute all of these.
#' @param probs A numeric vector of quantiles to compute.
#' @param boot_idx Logical indicating whether to use bootstrap samples in 
#' marginalization when computing quantiles. Default is `FALSE`. Currently 
#' only used in case where `condition` is `NULL` and 
#' @param ... Ignored.
#' @rdname predict
#' @export
predict.nhmm <- function(
    object, newdata, newdata2 = NULL, condition = NULL, 
    type = c("state", "response", "transition", "emission"),
    probs = c(0.025, 0.975), boot_idx = FALSE, ...) {
  # avoid CRAN check warnings due to NSE
  cols <- NULL
  type <- try(match.arg(type, several.ok = TRUE), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be {.val state}, {.val response}, 
    {.val transition}, {.val emission}, or a combination of these."
  )
  
  time <- object$time_variable
  id <- object$id_variable
  stopifnot_(
    !missing(newdata) && is.data.frame(newdata), 
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
  newdata <- copy(newdata)
  if (!is.null(newdata2)) {
    stopifnot_(
      is.data.frame(newdata2), 
      "Argument {.arg newdata2} must be a {.cls data.frame} object."
    )
    stopifnot_(
      identical(dim(newdata), dim(newdata2)), 
      "Data frames {.arg newdata} and {.arg newdata2} must have the same dimensions."
    )
    stopifnot_(
      identical(names(newdata), names(newdata2)), 
      "Data frames {.arg newdata} and {.arg newdata2} must have the same column names."
    )
    stopifnot_(
      identical(newdata[[id]], newdata2[[id]]),
      "Grouping variable {.var {id}} must be the same in both data frames."
    )
    stopifnot_(
      identical(newdata[[time]], newdata2[[time]]),
      "Time index variable {.var {time}} must be the same in both data frames."
    )
    newdata2 <- copy(newdata2)
  }
  
  responses <- object$responses
  S <- object$n_states
  M <- object$n_symbols
  C <- object$n_channels
  
  if (!is.null(condition)) {
    stopifnot_(
      all(condition %in% colnames(newdata)), 
      "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
    )
  }
  cond <- c(object$id_variable, condition)
  object <- update(object, newdata)
  if (inherits(object, "fanhmm")) {
    W <- update_W_for_fanhmm(object)
    object$W_A <- W$W_A
    object$W_B <- W$W_B
  }
  if (!is.null(newdata2)) {
    object2 <- update(object, newdata2)
    object2$boot <- NULL
    if (inherits(object2, "fanhmm")) {
      W <- update_W_for_fanhmm(object2)
      object2$W_A <- W$W_A
      object2$W_B <- W$W_B
    }
  }
  
  type_y <- "response" %in% type
  type_B <- "emission" %in% type
  type_z <- "state" %in% type
  type_A <- "transition" %in% type
  dy1 <- dy2 <- stats::setNames(vector("list", C), responses)
  dz1 <- dz2 <- NULL
  if (type_y || type_B) {
    for (y in responses) {
      idx <- rep(seq_row(object$data), each = S * M[y])
      dy1[[y]] <- setDT(object$data, key = c(id, time))[
        idx, cols, env = list(cols = as.list(cond))
      ]
      n <- nrow(dy1[[y]])
      state_names <- rep(object$state_names, times = M[y])
      symbol_names <- rep(object$symbol_names[[y]], each = S)
      set(dy1[[y]], j = "state", value = rep_len(state_names, n))
      set(dy1[[y]], j = y, value = rep_len(symbol_names, n))
      if (!is.null(newdata2)) {
        dy2[[y]] <- setDT(object2$data, key = c(id, time))[
          idx, cols, env = list(cols = as.list(cond))
        ]
        set(dy2[[y]], j = "state", value = rep_len(state_names, n))
        set(dy2[[y]], j = y, value = rep_len(symbol_names, n))
      }
    }
  }
  if (type_z || type_A) {
    idx <- rep(seq_row(object$data), each = S * S)
    dz1 <- setDT(object$data, key = c(id, time))[
      idx, cols, env = list(cols = as.list(cond))
    ]
    n <- nrow(dz1)
    state_names_from <- rep(object$state_names, times = S)
    state_names_to <- rep(object$state_names, each = S)
    set(dz1, j = "state_from", value = rep_len(state_names_from, n))
    set(dz1, j = "state_to", value = rep_len(state_names_to, n))
    if (!is.null(newdata2)) {
      dz2 <- setDT(object2$data, key = c(id, time))[
        idx, cols, env = list(cols = as.list(cond))
      ]
      set(dz2, j = "state_from", value = rep_len(state_names_from, n))
      set(dz2, j = "state_to", value = rep_len(state_names_to, n))
    }
  }
  
  obs <- create_obs(object)
  out <- .predict_nhmm(
    dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A
  )
  if (type_y) out_y <- out$y
  if (type_B) out_B <- out$B
  if (type_z) out_z <- out$z
  if (type_A) out_A <- out$A
  
  if (!is.null(newdata2)) {
    obs2 <- create_obs(object2)
    out <- .predict_nhmm(
      dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
    )
    if (type_y) prob_diff(out_y, out$y, responses)
    if (type_B) prob_diff(out_B, out$B, responses)
    if (type_z) prob_diff(out_z, out$z)
    if (type_A) prob_diff(out_A, out$A)
  }
  
  nsim <- length(object$boot$gamma_pi)
  if (nsim > 0 && length(probs) > 0) {
    boot_idx <- boot_idx && !is.null(object$boot$idx)
    if (type_y) {
      boot_y <- stats::setNames(vector("list", C), responses)
      for (y in responses) {
        boot_y[[y]] <- matrix(NA, nrow(out_y[[y]]), nsim)
      }
    }
    if (type_B) {
      boot_B <- stats::setNames(vector("list", C), responses)
      for (y in responses) {
        boot_B[[y]] <- matrix(NA, nrow(out_B[[y]]), nsim)
      }
    }
    if (type_z) boot_z <- matrix(NA, nrow(out_z), nsim)
    if (type_A) boot_A <- matrix(NA, nrow(out_A), nsim)
    for (i in seq_len(nsim)) {
      object$gammas$gamma_pi <- object$boot$gamma_pi[[i]]
      object$gammas$gamma_A <- object$boot$gamma_A[[i]]
      object$gammas$gamma_B <- object$boot$gamma_B[[i]]
      out <- .predict_nhmm(
        dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A
      )
      if (type_y) {
        for (y in responses) {
          boot_y[[y]][, i] <- out$y[[y]]$probability
        }
      }
      if (type_B) {
        for (y in responses) {
          boot_B[[y]][, i] <- out$B[[y]]$probability
        }
      }
      if (type_z) boot_z[, i] <- out$z$probability
      if (type_A) boot_A[, i] <- out$A$probability
      
      if (!is.null(newdata2)) {
        object2$gammas$gamma_pi <- object$boot$gamma_pi[[i]]
        object2$gammas$gamma_A <- object$boot$gamma_A[[i]]
        object2$gammas$gamma_B <- object$boot$gamma_B[[i]]
        out <- .predict_nhmm(
          dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
        )
        if (type_y) {
          for (y in responses) {
            boot_y[[y]][, i] <- boot_y[[y]][, i] - out$y[[y]]$probability
          }
        }
        if (type_B) {
          for (y in responses) {
            boot_B[[y]][, i] <- boot_B[[y]][, i] - out$B[[y]]$probability
          }
        }
        if (type_z) boot_z[, i] <- boot_z[, i] - out$z$probability
        if (type_A) boot_A[, i] <- boot_A[, i] - out$A$probability
      }
    }
    p <- paste0("q", 100 * probs)
    if (type_y) {
      for (y in responses) {
        q <- fast_quantiles(boot_y[[y]], probs)
        out_y[[y]][, (p) := data.table(q)]
      }
    }
    if (type_B) {
      for (y in responses) {
        q <- fast_quantiles(boot_B[[y]], probs)
        out_B[[y]][, (p) := data.table(q)]
      }
    }
    if (type_z) {
      q <- fast_quantiles(boot_z, probs)
      out_z[, (p) := data.table(q)]
    }
    if (type_A) {
      q <- fast_quantiles(boot_A, probs)
      out_A[, (p) := data.table(q)]
    }
  }
  out <- list(
    response = out_y, emission = out_B, state = out_z, transition = out_A
  )
  out[lengths(out) > 0]
}
#' @rdname predict
#' @export
predict.mnhmm <- function(
    object, newdata, newdata2 = NULL, condition = NULL, 
    type = c("state", "response", "transition", "emission"),
    probs = c(0.025, 0.975), boot_idx = FALSE, ...) {
  cols <- NULL
  type <- try(match.arg(type, several.ok = TRUE), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be {.val state}, {.val response}, 
    {.val transition}, {.val emission}, or a combination of these."
  )
  
  time <- object$time_variable
  id <- object$id_variable
  stopifnot_(
    !missing(newdata) && is.data.frame(newdata), 
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
  newdata <- copy(newdata)
  if (!is.null(newdata2)) {
    stopifnot_(
      is.data.frame(newdata2), 
      "Argument {.arg newdata2} must be a {.cls data.frame} object."
    )
    stopifnot_(
      identical(dim(newdata), dim(newdata2)), 
      "Data frames {.arg newdata} and {.arg newdata2} must have the same dimensions."
    )
    stopifnot_(
      identical(names(newdata), names(newdata2)), 
      "Data frames {.arg newdata} and {.arg newdata2} must have the same column names."
    )
    stopifnot_(
      identical(newdata[[id]], newdata2[[id]]),
      "Grouping variable {.var {id}} must be the same in both data frames."
    )
    stopifnot_(
      identical(newdata[[time]], newdata2[[time]]),
      "Time index variable {.var {time}} must be the same in both data frames."
    )
    newdata2 <- copy(newdata2)
  }
  
  responses <- object$responses
  S <- object$n_states
  M <- object$n_symbols
  C <- object$n_channels
  D <- object$n_clusters
  
  if (!is.null(condition)) {
    stopifnot_(
      all(condition %in% colnames(newdata)), 
      "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
    )
  }
  cond <- c(object$id_variable, condition)
  object <- update(object, newdata)
  if (inherits(object, "fanhmm")) {
    W <- update_W_for_fanhmm(object)
    object$W_A <- W$W_A
    object$W_B <- W$W_B
  }
  if (!is.null(newdata2)) {
    object2 <- update(object, newdata2)
    object2$boot <- NULL
    if (inherits(object2, "fanhmm")) {
      W <- update_W_for_fanhmm(object2)
      object2$W_A <- W$W_A
      object2$W_B <- W$W_B
    }
  }
  
  type_y <- "response" %in% type
  type_B <- "emission" %in% type
  type_z <- "state" %in% type
  type_A <- "transition" %in% type
  dy1 <- dy2 <- stats::setNames(vector("list", C), responses)
  dz1 <- dz2 <- NULL
  states <- object$state_names[[1]]
  if (type_y || type_B) {
    for (y in responses) {
      idx <- rep(seq_row(object$data), each = S * M[y])
      dy1[[y]] <- setDT(object$data, key = c(id, time))[
        idx, cols, env = list(cols = as.list(cond))
      ]
      n <- nrow(dy1[[y]])
      state_names <- rep(states, times = M[y])
      symbol_names <- rep(object$symbol_names[[y]], each = S)
      set(dy1[[y]], j = "state", value = rep_len(state_names, n))
      set(dy1[[y]], j = y, value = rep_len(symbol_names, n))
      if (!is.null(newdata2)) {
        dy2[[y]] <- setDT(object2$data, key = c(id, time))[
          idx, cols, env = list(cols = as.list(cond))
        ]
        set(dy2[[y]], j = "state", value = rep_len(state_names, n))
        set(dy2[[y]], j = y, value = rep_len(symbol_names, n))
      }
    }
  }
  if (type_z || type_A) {
    idx <- rep(seq_row(object$data), each = S * S)
    dz1 <- setDT(object$data, key = c(id, time))[
      idx, cols, env = list(cols = as.list(cond))
    ]
    n <- nrow(dz1)
    state_names_from <- rep(states, times = S)
    state_names_to <- rep(states, each = S)
    set(dz1, j = "state_from", value = rep_len(state_names_from, n))
    set(dz1, j = "state_to", value = rep_len(state_names_to, n))
    if (!is.null(newdata2)) {
      dz2 <- setDT(object2$data, key = c(id, time))[
        idx, cols, env = list(cols = as.list(cond))
      ]
      set(dz2, j = "state_from", value = rep_len(state_names_from, n))
      set(dz2, j = "state_to", value = rep_len(state_names_to, n))
    }
  }
  
  obs <- create_obs(object)
  out <- .predict_mnhmm(
    dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A
  )
  if (type_y) out_y <- out$y
  if (type_B) out_B <- out$B
  if (type_z) out_z <- out$z
  if (type_A) out_A <- out$A
  
  if (!is.null(newdata2)) {
    obs2 <- create_obs(object2)
    out <- .predict_mnhmm(
      dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
    )
    if (type_y) prob_diff(out_y, out$y, responses, clusters)
    if (type_B) prob_diff(out_B, out$B, responses, clusters)
    if (type_z) prob_diff(out_z, out$z, clusters)
    if (type_A) prob_diff(out_A, out$A, clusters)
  }
  
  nsim <- length(object$boot$gamma_pi)
  if (nsim > 0 && length(probs) > 0) {
    boot_idx <- boot_idx && !is.null(object$boot$idx)
    if (type_y) {
      boot_y <- stats::setNames(vector("list", C), responses)
      for (y in responses) {
        boot_y[[y]] <- array(NA, c(nrow(out_y[[1]][[y]]), D, nsim))
      }
    }
    if (type_B) {
      boot_B <- stats::setNames(vector("list", C), responses)
      for (y in responses) {
        boot_B[[y]] <- array(NA, c(nrow(out_B[[1]][[y]]), D, nsim))
      }
    }
    if (type_z) boot_z <- array(NA, c(nrow(out_z[[1]]), D, nsim))
    if (type_A) boot_A <- array(NA, c(nrow(out_A[[1]]), D, nsim))
    for (i in seq_len(nsim)) {
      object$gammas$gamma_pi <- object$boot$gamma_pi[[i]]
      object$gammas$gamma_A <- object$boot$gamma_A[[i]]
      object$gammas$gamma_B <- object$boot$gamma_B[[i]]
      object$gammas$gamma_omega <- object$boot$gamma_omega[[i]]
      out <- .predict_mnhmm(
        dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A
      )
      for (d in seq_len(D)) {
        if (type_y) {
          for (y in responses) {
            boot_y[[y]][, d, i] <- out$y[[d]][[y]]$probability
          }
        }
        if (type_B) {
          for (y in responses) {
            boot_B[[y]][, d, i] <- out$B[[d]][[y]]$probability
          }
        }
        if (type_z) boot_z[, d, i] <- out$z[[d]]$probability
        if (type_A) boot_A[, d, i] <- out$A[[d]]$probability
      }
      if (!is.null(newdata2)) {
        object2$gammas$gamma_pi <- object$boot$gamma_pi[[i]]
        object2$gammas$gamma_A <- object$boot$gamma_A[[i]]
        object2$gammas$gamma_B <- object$boot$gamma_B[[i]]
        object2$gammas$gamma_omega <- object$boot$gamma_omega[[i]]
        out <- .predict_mnhmm(
          dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
        )
        for (d in seq_len(D)) {
          if (type_y) {
            for (y in responses) {
              boot_y[[y]][, d, i] <- boot_y[[y]][, d, i] - out$y[[d]][[y]]$probability
            }
          }
          if (type_B) {
            for (y in responses) {
              boot_B[[y]][, d, i] <- boot_B[[y]][, d, i] - out$B[[d]][[y]]$probability
            }
          }
          if (type_z) boot_z[, d, i] <- boot_z[, d, i] - out$z[[d]]$probability
          if (type_A) boot_A[, d, i] <- boot_A[, d, i] - out$A[[d]]$probability
        }
      }
    }
    p <- paste0("q", 100 * probs)
    for (d in seq_len(D)) {
      if (type_y) {
        for (y in responses) {
          q <- fast_quantiles(boot_y[[y]][, d, ], probs)
          out_y[[d]][[y]][, (p) := data.table(q)]
        }
      }
      if (type_B) {
        for (y in responses) {
          q <- fast_quantiles(boot_B[[y]][, d, ], probs)
          out_B[[d]][[y]][, (p) := data.table(q)]
        }
      }
      if (type_z) {
        q <- fast_quantiles(boot_z[, d, ], probs)
        out_z[[d]][, (p) := data.table(q)]
      }
      if (type_A) {
        q <- fast_quantiles(boot_A[, d, ], probs)
        out_A[[d]][, (p) := data.table(q)]
      }
    }
  }
  out <- list(
    response = out_y, emission = out_B, state = out_z, transition = out_A
  )
  out[lengths(out) > 0]
}
