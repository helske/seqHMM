#' Compute marginal and conditional probabilities from joint distributions 
#' obtained from predict based on forward predictions
#' @noRd
predict_y_B <- function(x, newdata, response, cond, only_estimates = TRUE, 
                        cluster = NULL) {
  # avoid CRAN check warnings due to NSE
  estimate <- probability <- NULL
  set(newdata, j = "estimate", value = x)
  cond_all <- c(cond, cluster, "state", response)
  # Compute joint distribution P(Y_t, Z_t)
  d <- newdata[, list(probability = mean(estimate)), by = cond_all]
  setorderv(d, cond_all)
  if (only_estimates) {
    out_y <- d[, list(probability = sum(probability)), 
               by = c(cluster, response, cond)]$probability
    out_B <- d[, "probability" := list(probability / sum(probability)), 
               by = c(cluster, "state", cond)]$probability
  } else {
    out_y <- d[, list(probability = sum(probability)), 
               by = c(cluster, response, cond)]
    out_B <- d[, "probability" := list(probability / sum(probability)), 
               by = c(cluster, "state", cond)]
  }
  list(y = out_y, B = out_B)
}
predict_z_A <- function(x, newdata, cond, only_estimates = TRUE, cluster = NULL) {
  # avoid CRAN check warnings due to NSE
  estimate <- probability <- NULL
  set(newdata, j = "estimate", value = x)
  cond_all <- c(cond, cluster, "state_from", "state_to")
  # Compute joint distribution P(Z_t-1, Z_t)
  d <- newdata[, list(probability = mean(estimate)), by = cond_all]
  setorderv(d, cond_all)
  if (only_estimates) {
    out_z <- d[, list(probability = sum(probability)), 
               by = c(cluster, "state_to", cond)]$probability
    out_A <- d[, "probability" := list(probability / sum(probability)), 
               by = c(cluster, "state_from", cond)]$probability
  } else {
    out_z <- d[, list(probability = sum(probability)), by = c(cluster, "state_to", cond)]
    out_A <- d[, "probability" := list(probability / sum(probability)), 
               by = c(cluster, "state_from", cond)]
  }
  list(z = out_z, A = out_A)
}

.predict <- function(
    dy, dz, condition, obs, object, type_y, type_B, type_z, type_A, 
    only_estimates = TRUE) {
  
  if (inherits(object, "fanhmm")) {
    if (inherits(object, "nhmm")) {
      out <- Rcpp_predict_fanhmm(
        obs, object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, 
        io(object$X_pi), io(object$X_A), io(object$X_B), 
        iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
        object$gammas$gamma_pi, object$gammas$gamma_A, object$gammas$gamma_B,
        object$prior_y0, object$W_X_B, object$W_A, object$W_B
      )
    } else {
      out <- Rcpp_predict_mfanhmm(
        obs, object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, object$X_omega,
        io(object$X_pi), io(object$X_A), io(object$X_B), io(object$X_omega),
        iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
        object$gammas$gamma_pi, object$gammas$gamma_A, object$gammas$gamma_B,
        object$gammas$gamma_omega,
        object$prior_y0, object$W_X_B, object$W_A, object$W_B
      )
    }
  } else {
    if (inherits(object, "nhmm")) {
      out <- Rcpp_predict_nhmm(
        obs, object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, 
        io(object$X_pi), io(object$X_A), io(object$X_B), 
        iv(object$X_A), iv(object$X_B), tv(object$X_A), tv(object$X_B),
        object$gammas$gamma_pi, object$gammas$gamma_A, object$gammas$gamma_B
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
  }
  results_y <- results_B <- stats::setNames(
    vector("list", object$n_channels), object$responses
  )
  results_z <- results_A <- NULL
  if (inherits(object, "mnhmm")) {
    cluster <- "cluster"
  } else {
    cluster <- NULL
  }
  if (type_y || type_B) {
    for (y in seq_along(object$responses)) {
      d <- predict_y_B(
        unlist(out$observations[y, ,]), dy[[y]], object$responses[y], 
        condition, only_estimates, cluster
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
    d <- predict_z_A(unlist(out$states), dz, condition, only_estimates, cluster)
    if (type_z) {
      results_z <- d$z
    }
    if (type_A) {
      results_A <- d$A
    }
  }
  list(y = results_y, B = results_B, z = results_z, A = results_A)
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
      idx <- rep(seq_len(nrow(object$data)), each = S * M[y])
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
    idx <- rep(seq_len(nrow(object$data)), each = S * S)
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
  
  obs <- create_obsArray(object)
  out <- .predict(
    dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A,
    only_estimates = FALSE
  )
  if (type_y) out_y <- out$y
  if (type_B) out_B <- out$B
  if (type_z) out_z <- out$z
  if (type_A) out_A <- out$A
  
  if (!is.null(newdata2)) {
    obs2 <- create_obsArray(object2)
    out <- .predict(
      dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
    )
    if (type_y) {
      for (y in responses) {
        out_y[[y]][, probability := probability - out$y[[y]]]
      }
    } 
    if (type_B) {
      for (y in responses) {
        out_B[[y]][, probability := probability - out$B[[y]]]
      }
    }
    if (type_z) out_z[, probability := probability - out$z]
    if (type_A) out_A[, probability := probability - out$A]
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
      out <- .predict(
        dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A
      )
      if (type_y) {
        for (y in responses) {
          boot_y[[y]][, i] <- out$y[[y]]
        }
      }
      if (type_B) {
        for (y in responses) {
          boot_B[[y]][, i] <- out$B[[y]]
        }
      }
      if (type_z) boot_z[, i] <- out$z
      if (type_A) boot_A[, i] <- out$A
      
      if (!is.null(newdata2)) {
        object2$gammas$gamma_pi <- object$boot$gamma_pi[[i]]
        object2$gammas$gamma_A <- object$boot$gamma_A[[i]]
        object2$gammas$gamma_B <- object$boot$gamma_B[[i]]
        out <- .predict(
          dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
        )
        if (type_y) {
          for (y in responses) {
            boot_y[[y]][, i] <- boot_y[[y]][, i] - out$y[[y]]
          }
        }
        if (type_B) {
          for (y in responses) {
            boot_B[[y]][, i] <- boot_B[[y]][, i] - out$B[[y]]
          }
        }
        if (type_z) boot_z[, i] <- boot_z[, i] - out$z
        if (type_A) boot_A[, i] <- boot_A[, i] - out$A
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
  if (type_y || type_B) {
    for (y in responses) {
      idx <- rep(seq_len(nrow(newdata)), each = D * S * M[y])
      dy1[[y]] <- setDT(newdata, key = c(id, time))[
        idx, cols, env = list(cols = as.list(cond))
      ]
      n <- nrow(dy1[[y]])
      cluster_names <- rep(rep(object$cluster_names, each = S), times = M[y])
      state_names <- rep(unlist(object$state_names), times = M[y])
      symbol_names <- rep(object$symbol_names[[y]], each = D * S)
      set(dy1[[y]], j = "cluster", value = rep_len(cluster_names, n))
      set(dy1[[y]], j = "state", value = rep_len(state_names, n))
      set(dy1[[y]], j = y, value = rep_len(symbol_names, n))
      if (!is.null(newdata2)) {
        dy2[[y]] <- setDT(newdata2, key = c(id, time))[
          idx, cols, env = list(cols = as.list(cond))
        ]
        set(dy2[[y]], j = "cluster", value = rep_len(cluster_names, n))
        set(dy2[[y]], j = "state", value = rep_len(state_names, n))
        set(dy2[[y]], j = y, value = rep_len(symbol_names, n))
      }
    }
  }
  if (type_z || type_A) {
    idx <- rep(seq_len(nrow(newdata)), each = D * S * S)
    dz1 <- setDT(newdata, key = c(id, time))[
      idx, cols, env = list(cols = as.list(cond))
    ]
    n <- nrow(dz1)
    cluster_names <- rep(rep(object$cluster_names, each = S), times = S)
    state_names_from <- rep(unlist(object$state_names), times = S)
    state_names_to <- rep(unlist(object$state_names), each = D * S)
    set(dz1, j = "cluster", value = rep_len(cluster_names, n))
    set(dz1, j = "state_from", value = rep_len(state_names_from, n))
    set(dz1, j = "state_to", value = rep_len(state_names_to, n))
    if (!is.null(newdata2)) {
      dz2 <- setDT(newdata2, key = c(id, time))[
        idx, cols, env = list(cols = as.list(cond))
      ]
      set(dz2, j = "cluster", value = rep_len(cluster_names, n))
      set(dz2, j = "state_from", value = rep_len(state_names_from, n))
      set(dz2, j = "state_to", value = rep_len(state_names_to, n))
    }
  }
  
  obs <- create_obsArray(object)
  out <- .predict(
    dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A,
    only_estimates = FALSE
  )
  if (type_y) out_y <- out$y
  if (type_B) out_B <- out$B
  if (type_z) out_z <- out$z
  if (type_A) out_A <- out$A
  
  if (!is.null(newdata2)) {
    obs2 <- create_obsArray(object2)
    out <- .predict(
      dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
    )
    if (type_y) {
      for (y in responses) {
        out_y[[y]][, probability := probability - out$y[[y]]]
      }
    } 
    if (type_B) {
      for (y in responses) {
        out_B[[y]][, probability := probability - out$B[[y]]]
      }
    }
    if (type_z) out_z[, probability := probability - out$z]
    if (type_A) out_A[, probability := probability - out$A]
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
      object$gammas$gamma_omega <- object$boot$gamma_omega[[i]]
      out <- .predict(
        dy1, dz1, condition, obs, object, type_y, type_B, type_z, type_A
      )
      if (type_y) {
        for (y in responses) {
          boot_y[[y]][, i] <- out$y[[y]]
        }
      }
      if (type_B) {
        for (y in responses) {
          boot_B[[y]][, i] <- out$B[[y]]
        }
      }
      if (type_z) boot_z[, i] <- out$z
      if (type_A) boot_A[, i] <- out$A
      
      if (!is.null(newdata2)) {
        object2$gammas$gamma_pi <- object$boot$gamma_pi[[i]]
        object2$gammas$gamma_A <- object$boot$gamma_A[[i]]
        object2$gammas$gamma_B <- object$boot$gamma_B[[i]]
        object2$gammas$gamma_omega <- object$boot$gamma_omega[[i]]
        out <- .predict(
          dy2, dz2, condition, obs2, object2, type_y, type_B, type_z, type_A
        )
        if (type_y) {
          for (y in responses) {
            boot_y[[y]][, i] <- boot_y[[y]][, i] - out$y[[y]]
          }
        }
        if (type_B) {
          for (y in responses) {
            boot_B[[y]][, i] <- boot_B[[y]][, i] - out$B[[y]]
          }
        }
        if (type_z) boot_z[, i] <- boot_z[, i] - out$z
        if (type_A) boot_A[, i] <- boot_A[, i] - out$A
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
