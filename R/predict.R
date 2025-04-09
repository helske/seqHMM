#' Compute marginal and conditional probabilities from joint distributions obtained from predict
#' @noRd
compute_joint <- function(x, newdata, type, cond, include) {
  # avoid CRAN check warnings due to NSE
  estimate <- probability <- NULL
  cond_obs <- cond$obs
  cond_state <- cond$state
  cond_both <- cond$both
  set(newdata, j = "estimate", value = stats::na.omit(c(x)))
  # Compute joint distribution P(Y, Z)
  if (!is.null(include)) {
    newdata <- newdata[include]
  }
  d <- newdata[, list(probability = mean(estimate)), by = cond_both]
  d_obs <- d_state <- d_cond <- NULL
  
  # P(Y)
  if ("observations" %in% type) {
    d_obs <- d[, list(probability = sum(probability)), by = cond_obs]
    setorderv(d_obs, cond_obs)
  }
  # P(Z)
  if ("states" %in% type) {
    d_state <- d[, list(probability = sum(probability)), by = cond_state]
    setorderv(d_state, cond_state)
  }
  # P(Y | Z)
  if ("conditionals" %in% type) {
    d_cond <- d[, "probability" := list(probability / sum(probability)), by = cond_state]
    setorderv(d_cond, cond_state)
  }
  list(observations = d_obs, states = d_state, conditionals = d_cond)
}
#' Predictions from Non-homogeneous Hidden Markov Models
#'
#' @rdname predict
#' @param object An object of class `nhmm`, `mnhmm`, or `fanhmm`.
#' @param newdata A data frame used for computing the predictions.
#' @param newdata2 An optional data frame for predictions, in which case the 
#' estimates are differences between predictions using `newdata` and `newdata2`.
#' @param condition An optional vector of variable names used for conditional 
#' predictions.
#' @param type A character string specifying the type of predictions. Possible 
#' values are `"observations"` (marginal probabilities of observations), 
#' `"states"` (marginals of states), and `"conditionals"` (emission probabilities).
#' Default computes all of these.
#' @param probs A numeric vector of quantiles to compute.
#' @param boot_idx Logical indicating whether to use bootstrap samples in 
#' marginalization when computing quantiles. Default is `FALSE`.
#' @param ... Ignored.
#' @export
predict.nhmm <- function(
    object, newdata, newdata2 = NULL, condition = NULL, 
    type = c("observations", "states", "conditionals"),
    probs = c(0.025, 0.975), boot_idx = FALSE, ...) {
  cols <- NULL
  type <- try(match.arg(type, several.ok = TRUE), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be {.val observations},
    {.val states}, {.val conditionals}, or a combination of these."
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
  
  state_names <- object$state_names
  symbol_names <- object$symbol_names
  responses <- object$responses
  S <- object$n_states
  M <- object$n_symbols
  
  if (!is.null(condition)) {
    stopifnot_(
      all(condition %in% colnames(newdata)), 
      "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
    )
  }
  conditions <- list(
    base = condition,
    obs = c(responses, condition),
    state = c("state", condition),
    both = c("state", responses, condition)
  )
  if (is.null(conditions$base)) conditions$base <- id
  object <- update(object, newdata)
  obs <- create_obsArray(object)
  state_factor <- factor(rep(state_names, times = M), levels = state_names)
  obs_factor <- factor(rep(symbol_names, each = S), levels = symbol_names)
  newdata <- setDT(newdata, key = c(id, time))[
    rep(seq_len(nrow(newdata)), each = S * M), 
    cols, env = list(cols = as.list(conditions$base))
  ]
  set(newdata, j = "state", value = rep(state_factor, length.out = nrow(newdata)))
  set(newdata, j = responses, value = rep(obs_factor, length.out = nrow(newdata)))
  include <- NULL
  
  out <- unlist(predict_nhmm(
    obs, object$sequence_lengths, object$n_symbols, 
    object$X_pi, object$X_A, object$X_B, 
    io(object$X_pi), io(object$X_A), io(object$X_B), 
    iv(object$X_A), iv(object$X_B),
    tv(object$X_A), tv(object$X_B),
    object$etas$pi, object$etas$A, object$etas$B
  ))
  d_mean <- compute_joint(out, newdata, type, conditions, include)
  
  if (!is.null(newdata2)) {
    object2 <- update(object, newdata2)
    newdata2 <- setDT(newdata2, key = c(id, time))[
      rep(seq_len(nrow(newdata2)), each = S * M), 
      cols, env = list(cols = as.list(conditions$base))
    ]
    set(newdata2, j = "state", value = rep(state_factor, length.out = nrow(newdata2)))
    set(newdata2, j = responses, value = rep(obs_factor, length.out = nrow(newdata2)))
    obs2 <- create_obsArray(object2)
    out2 <- unlist(predict_nhmm( 
      obs2, object2$sequence_lengths, object2$n_symbols, 
      object2$X_pi, object2$X_A, object2$X_B, 
      io(object2$X_pi), io(object2$X_A), io(object2$X_B), 
      iv(object2$X_A), iv(object2$X_B),
      tv(object2$X_A), tv(object2$X_B),
      object2$etas$pi, object2$etas$A, object2$etas$B
    ))
    d <- compute_joint(out2, newdata2, type, conditions, include)
    for(i in names(d_mean)) {
      if(!is.null(d_mean[[i]])) {
        d_mean[[i]]$probability <- d_mean[[i]]$probability - d[[i]]$probability
      }
    }
  }
  
  nsim <- length(object$boot$gamma_pi)
  if (nsim > 0 && length(probs) > 0) {
    d_boot <- vector("list", nsim)
    boot_idx <- boot_idx && !is.null(object$boot$idx)
    for (i in seq_len(nsim)) {
      out <- simplify2array(boot_predict_nhmm( 
        obs, object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, 
        io(object$X_pi), io(object$X_A), io(object$X_B), 
        iv(object$X_A), iv(object$X_B),
        tv(object$X_A), tv(object$X_B),
        object$etas$pi, object$etas$A, object$etas$B,
        object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]]
      ))
      if (boot_idx) {
        d_boot[[i]] <- compute_joint(
          out[, , , object$boot$idx[, i], drop = FALSE], 
          newdata, type, conditions, include
        )
      } else {
        d_boot[[i]] <- compute_joint(out, newdata, type, conditions, include)
      }
      if (!is.null(newdata2)) {
        out2 <- simplify2array(boot_predict_nhmm( 
          obs2, object2$sequence_lengths, object2$n_symbols, 
          object2$X_pi, object2$X_A, object2$X_B, 
          io(object2$X_pi), io(object2$X_A), io(object2$X_B), 
          iv(object2$X_A), iv(object2$X_B),
          tv(object2$X_A), tv(object2$X_B),
          object2$etas$pi, object2$etas$A, object2$etas$B,
          object2$boot$gamma_pi[[i]], object2$boot$gamma_A[[i]], object2$boot$gamma_B[[i]]
        ))
        if (boot_idx) {
          d <- compute_joint(
            out2[, , , object$boot$idx[, i], drop = FALSE], 
            newdata2, type, conditions, include
          )
        } else {
          d <- compute_joint(
            out2, newdata2, type, conditions, include)
        }
        for(j in names(d_mean)) {
          if(!is.null(d_mean[[j]])) {
            d_boot[[i]][[j]]$probability <- d_boot[[i]][[j]]$probability - d[[j]]$probability
          }
        }
      }
    }
    for(i in names(d_mean)) {
      if(!is.null(d_mean[[i]])) {
        di <- do.call(rbind, lapply(d_boot, \(x) x[[i]]$probability))
        qs <- t(apply(di, 2, stats::quantile, probs = probs))
        colnames(qs) <- paste0("q", 100 * probs)
        d_mean[[i]] <- cbind(d_mean[[i]], qs)
      }
    }
  }
  lapply(d_mean[lengths(d_mean) > 0], as.data.frame)
}
#' @rdname predict
#' @export
predict.fanhmm <- function(
    object, newdata, newdata2 = NULL, condition = NULL, 
    type = c("observations", "states", "conditionals"),
    probs = c(0.025, 0.975), boot_idx = FALSE, ...) {
  cols <- NULL
  stopifnot_(
    object$n_channels == 1L,
    "Multichannel FAN-HMM is not yet supported."
  )
  type <- try(match.arg(type, several.ok = TRUE), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be {.val observations},
    {.val states}, {.val conditionals}, or a combination of these."
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
  
  state_names <- object$state_names
  symbol_names <- object$symbol_names
  responses <- object$responses
  S <- object$n_states
  M <- object$n_symbols
  
  if (!is.null(condition)) {
    stopifnot_(
      all(condition %in% colnames(newdata)), 
      "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
    )
  }
  conditions <- list(
    base = condition,
    obs = c(responses, condition),
    state = c("state", condition),
    both = c("state", responses, condition)
  )
  if (is.null(conditions$base)) conditions$base <- id
  object <- update(object, newdata)
  obs <- create_obsArray(object)
  state_factor <- factor(rep(state_names, times = M), levels = state_names)
  obs_factor <- factor(rep(symbol_names, each = S), levels = symbol_names)
  newdata <- setDT(newdata, key = c(id, time))[
    rep(seq_len(nrow(newdata)), each = S * M), 
    cols, env = list(cols = as.list(c(time, conditions$base)))
  ]
  set(newdata, j = "state", value = rep(state_factor, length.out = nrow(newdata)))
  set(newdata, j = responses, value = rep(obs_factor, length.out = nrow(newdata)))
  if (!is.null(object$autoregression_formula)) {
    include <- which(newdata[[time]] != min(newdata[[time]]))
  } else {
    include <- NULL
  }
  newdata[, time := NULL, env = list(time = time)]
  
  W <- update_W_for_fanhmm(object)
  out <- unlist(predict_fanhmm( 
    obs, object$sequence_lengths, object$n_symbols, 
    object$X_pi, object$X_A, object$X_B, 
    io(object$X_pi), io(object$X_A), io(object$X_B), 
    iv(object$X_A), iv(object$X_B),
    tv(object$X_A), tv(object$X_B),
    object$etas$pi, object$etas$A, object$etas$B,
    W$W_A, W$W_B, !is.null(object$autoregression_formula)
  ))
  d_mean <- compute_joint(out, newdata, type, conditions, include)
  if (!is.null(newdata2)) {
    object2 <- update(object, newdata2)
    newdata2 <- setDT(newdata2, key = c(id, time))[
      rep(seq_len(nrow(newdata2)), each = S * M), 
      cols, env = list(cols = as.list(conditions$base))
    ]
    set(newdata2, j = "state", value = rep(state_factor, length.out = nrow(newdata2)))
    set(newdata2, j = responses, value = rep(obs_factor, length.out = nrow(newdata2)))
    obs2 <- create_obsArray(object2)
    W2 <- update_W_for_fanhmm(object2)
    out2 <- unlist(predict_fanhmm( 
      obs2, object2$sequence_lengths, object2$n_symbols, 
      object2$X_pi, object2$X_A, object2$X_B, 
      io(object2$X_pi), io(object2$X_A), io(object2$X_B), 
      iv(object2$X_A), iv(object2$X_B),
      tv(object2$X_A), tv(object2$X_B),
      object2$etas$pi, object2$etas$A, object2$etas$B,
      W2$W_A, W2$W_B, !is.null(object2$autoregression_formula)
    ))
    d <- compute_joint(out2, newdata2, type, conditions, include)
    for(i in names(d_mean)) {
      if(!is.null(d_mean[[i]])) {
        d_mean[[i]]$probability <- d_mean[[i]]$probability - d[[i]]$probability
      }
    }
  }
  
  nsim <- length(object$boot$gamma_pi)
  if (nsim > 0 && length(probs) > 0) {
    d_boot <- vector("list", nsim)
    boot_idx <- boot_idx && !is.null(object$boot$idx)
    for (i in seq_len(nsim)) {
      out <- simplify2array(boot_predict_fanhmm( 
        obs, object$sequence_lengths, object$n_symbols, 
        object$X_pi, object$X_A, object$X_B, 
        io(object$X_pi), io(object$X_A), io(object$X_B), 
        iv(object$X_A), iv(object$X_B),
        tv(object$X_A), tv(object$X_B),
        object$etas$pi, object$etas$A, object$etas$B,
        W$W_A, W$W_B,
        object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]],
        !is.null(object$autoregression_formula)
      ))
      if (boot_idx) {
        d_boot[[i]] <- compute_joint(
          out[, , , object$boot$idx[, i], drop = FALSE], 
          newdata, type, conditions, include
        )
      } else {
        d_boot[[i]] <- compute_joint(
          out, newdata, type, conditions, include
        )
      }
      if (!is.null(newdata2)) {
        out2 <- simplify2array(boot_predict_fanhmm( 
          obs2, object2$sequence_lengths, object2$n_symbols, 
          object2$X_pi, object2$X_A, object2$X_B, 
          io(object2$X_pi), io(object2$X_A), io(object2$X_B), 
          iv(object2$X_A), iv(object2$X_B),
          tv(object2$X_A), tv(object2$X_B),
          object2$etas$pi, object2$etas$A, object2$etas$B,
          W2$W_A, W2$W_B,
          object2$boot$gamma_pi[[i]], object2$boot$gamma_A[[i]], object2$boot$gamma_B[[i]],
          !is.null(object2$autoregression_formula)
        ))
        if (boot_idx) {
          d <- compute_joint(
            out2[, , , object$boot$idx[, i], drop = FALSE], 
            newdata2, type, conditions, include
          )
        } else {
          d <- compute_joint(out2, newdata2, type, conditions, include)
        }
        for(j in names(d_mean)) {
          if(!is.null(d_mean[[j]])) {
            d_boot[[i]][[j]]$probability <- d_boot[[i]][[j]]$probability - d[[j]]$probability
          }
        }
      }
    }
    for(i in names(d_mean)) {
      if(!is.null(d_mean[[i]])) {
        di <- do.call(rbind, lapply(d_boot, \(x) x[[i]]$probability))
        qs <- t(apply(di, 2, quantileq, probs = probs))
        d_mean[[i]] <- cbind(d_mean[[i]], qs)
      }
    }
  }
  lapply(d_mean[lengths(d_mean) > 0], as.data.table)
}
