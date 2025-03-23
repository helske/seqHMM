#' Compute marginal and conditional probabilities from joint distributions obtained from predict
#' @noRd
compute_joint <- function(x, newdata, type, cond, state_names, symbol_names, responses) {
  # avoid CRAN check warnings due to NSE
  ..cond_base <- estimate <- probability <- NULL
  S <- length(state_names)
  M <- length(symbol_names)
  cond_obs <- cond$obs
  cond_state <- cond$state
  cond_both <- cond$both
  cond_base <- cond$base
  # Convert inputs to factors
  state_factor <- factor(rep(state_names, times = M), levels = state_names)
  obs_factor <- factor(rep(symbol_names, each = S), levels = symbol_names)
  newdata <- setDT(newdata)[rep(seq_len(nrow(newdata)), each = S * M), 
                            ..cond_base]
  set(newdata, j = "state", value = rep(state_factor, length.out = nrow(newdata)))
  set(newdata, j = responses, value = rep(obs_factor, length.out = nrow(newdata)))
  set(newdata, j = "estimate", value = stats::na.omit(c(x)))
  # Compute joint distribution P(Y, Z)
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
    probs = c(0.025, 0.975), 
    boot_idx = FALSE, ...) {
  
  type <- match.arg(type, several.ok = TRUE)
  
  time <- object$time_variable
  id <- object$id_variable
  stopifnot_(
    !missing(newdata) & is.data.frame(newdata), 
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
  obsArray <- create_obsArray(object)
  stopifnot_(
    object$n_channels == 1L,
    "Predictions for multichannel NHMMs are not yet supported."
  )
  
  out <- unlist(predict_nhmm_singlechannel( 
    object$etas$pi, object$X_pi,
    object$etas$A, object$X_A, 
    object$etas$B, object$X_B, 
    obsArray,
    object$sequence_lengths, 
    attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
    attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
    attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv")
  ))
  d_mean <- compute_joint(out, newdata, type, conditions, state_names, symbol_names, responses)
  
  if (!is.null(newdata2)) {
    object2 <- update(object, newdata2)
    obs2 <- create_obsArray(object2)
    out2 <- unlist(predict_nhmm_singlechannel( 
      object2$etas$pi, object2$X_pi,
      object2$etas$A, object2$X_A, 
      object2$etas$B, object2$X_B, 
      array(obs2[1, , ], dim(obs2)[2:3]),
      object2$sequence_lengths, 
      attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
      attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
      attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv")
    ))
    d <- compute_joint(
      out2, newdata2, type, conditions, state_names, symbol_names, 
      responses
    )
    for(i in names(d_mean)) {
      if(!is.null(d_mean[[i]])) {
        d_mean[[i]]$probability <- d_mean[[i]]$probability - d[[i]]$probability
      }
    }
  }
  
  nsim <- length(object$boot$gamma_pi)
  if (nsim > 0 && length(probs) > 0) {
    d_boot <- vector("list", nsim)
    boot_idx <- boot_idx & !is.null(object$boot$idx)
    for (i in seq_len(nsim)) {
      out <- simplify2array(boot_predict_nhmm_singlechannel( 
        object$etas$pi, object$X_pi,
        object$etas$A, object$X_A, 
        object$etas$B, object$X_B, 
        obsArray,
        object$sequence_lengths, 
        attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
        attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
        attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
        object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]]
      ))
      if (boot_idx) {
        d_boot[[i]] <- compute_joint(
          out[, , , object$boot$idx[, i], drop = FALSE], 
          newdata, type, conditions, state_names, symbol_names, responses
        )
      } else {
        d_boot[[i]] <- compute_joint(
          out, newdata, type, conditions, state_names, symbol_names, 
          responses
        )
      }
      if (!is.null(newdata2)) {
        out2 <- simplify2array(boot_predict_nhmm_singlechannel( 
          object2$etas$pi, object2$X_pi,
          object2$etas$A, object2$X_A, 
          object2$etas$B, object2$X_B, 
          array(obs2[1, , ], dim(obs2)[2:3]),
          object2$sequence_lengths, 
          attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
          attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
          attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
          object2$boot$gamma_pi[[i]], object2$boot$gamma_A[[i]], object2$boot$gamma_B[[i]]
        ))
        if (boot_idx) {
          d <- compute_joint(
            out2[, , , object$boot$idx[, i], drop = FALSE], 
            newdata2, type, conditions, state_names, symbol_names,
            responses
          )
        } else {
          d <- compute_joint(
            out2, newdata2, type, conditions, state_names, symbol_names,
            responses
          )
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
        di <- do.call(rbind, lapply(d_boot, function(x) x[[i]]$probability))
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
    probs = c(0.025, 0.975), 
    boot_idx = FALSE, ...) {
  
  type <- match.arg(type, several.ok = TRUE)
  
  time <- object$time_variable
  id <- object$id_variable
  stopifnot_(
    !missing(newdata) & is.data.frame(newdata), 
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
  obs <- create_obsArray(object, FALSE) # don't set first obs to missing
  stopifnot_(
    object$n_channels == 1L,
    "Multichannel FAN-HMM is not yet supported."
  )
  W <- update_W_for_fanhmm(object, newdata)
  
  out <- unlist(predict_fanhmm_singlechannel( 
    object$etas$pi, object$X_pi,
    object$etas$A, object$X_A, 
    object$etas$B, object$X_B, 
    obs, object$sequence_lengths, 
    attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
    attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
    attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
    W$W_A, W$W_B, !is.null(object$autoregression_formula)
  ))
  d_mean <- compute_joint(out, newdata, type, conditions, state_names, symbol_names, responses)
  
  if (!is.null(newdata2)) {
    object2 <- update(object, newdata2)
    obs2 <- create_obsArray(object2, FALSE) # don't set first obs to missing
    W2 <- update_W_for_fanhmm(object2, newdata2)
    out2 <- unlist(predict_fanhmm_singlechannel( 
      object2$etas$pi, object2$X_pi,
      object2$etas$A, object2$X_A, 
      object2$etas$B, object2$X_B, 
      obs2, object2$sequence_lengths, 
      attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
      attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
      attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
      W2$W_A, W2$W_B, !is.null(object2$autoregression_formula)
    ))
    d <- compute_joint(
      out2, newdata2, type, conditions, state_names, symbol_names, 
      responses
    )
    for(i in names(d_mean)) {
      if(!is.null(d_mean[[i]])) {
        d_mean[[i]]$probability <- d_mean[[i]]$probability - d[[i]]$probability
      }
    }
  }
  
  nsim <- length(object$boot$gamma_pi)
  if (nsim > 0 && length(probs) > 0) {
    d_boot <- vector("list", nsim)
    boot_idx <- boot_idx & !is.null(object$boot$idx)
    for (i in seq_len(nsim)) {
      out <- simplify2array(boot_predict_fanhmm_singlechannel( 
        object$etas$pi, object$X_pi,
        object$etas$A, object$X_A, 
        object$etas$B, object$X_B, 
        obs, object$sequence_lengths, 
        attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
        attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
        attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
        W$W_A, W$W_B,
        object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]],
        !is.null(object$autoregression_formula)
      ))
      if (boot_idx) {
        d_boot[[i]] <- compute_joint(
          out[, , , object$boot$idx[, i], drop = FALSE], 
          newdata, type, conditions, state_names, symbol_names, responses
        )
      } else {
        d_boot[[i]] <- compute_joint(
          out, newdata, type, conditions, state_names, symbol_names, 
          responses
        )
      }
      if (!is.null(newdata2)) {
        out2 <- simplify2array(boot_predict_fanhmm_singlechannel( 
          object2$etas$pi, object2$X_pi,
          object2$etas$A, object2$X_A, 
          object2$etas$B, object2$X_B, 
          obs2, object2$sequence_lengths, 
          attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
          attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
          attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
          W2$W_A, W2$W_B,
          object2$boot$gamma_pi[[i]], object2$boot$gamma_A[[i]], object2$boot$gamma_B[[i]],
          !is.null(object2$autoregression_formula)
        ))
        if (boot_idx) {
          d <- compute_joint(
            out2[, , , object$boot$idx[, i], drop = FALSE], 
            newdata2, type, conditions, state_names, symbol_names,
            responses
          )
        } else {
          d <- compute_joint(
            out2, newdata2, type, conditions, state_names, symbol_names,
            responses
          )
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
        di <- do.call(rbind, lapply(d_boot, function(x) x[[i]]$probability))
        qs <- t(apply(di, 2, stats::quantile, probs = probs))
        colnames(qs) <- paste0("q", 100 * probs)
        d_mean[[i]] <- cbind(d_mean[[i]], qs)
      }
    }
  }
  lapply(d_mean[lengths(d_mean) > 0], as.data.frame)
}
