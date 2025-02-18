#' Compute marginal and conditional probabilities from joint distributions obtained from predict
#' @noRd 
#' @importFrom data.table set setDT := copy
compute_joint <- function(x, newdata, type, cond, state_names, symbol_names, response_names) {
  S <- length(state_names)
  M <- length(symbol_names)
  cond_obs <- cond$obs
  cond_state <- cond$state
  cond_both <- cond$both
  cond_base <- cond$base
  # Convert inputs to factors
  state_factor <- factor(rep(state_names, times = M), levels = state_names)
  obs_factor <- factor(rep(symbol_names, each = S), levels = symbol_names)
  newdata <- setDT(newdata)[rep(seq_len(nrow(newdata)), each = S * M), cond_base, with = FALSE, drop = FALSE]
  set(newdata, j = "state", value = rep(state_factor, length.out = nrow(newdata)))
  set(newdata, j = response_names, value = rep(obs_factor, length.out = nrow(newdata)))
  set(newdata, j = "estimate", value = na.omit(c(x)))
  # Compute joint distribution P(Y, Z)
  d <- newdata[, .(probability = mean(estimate)), by = cond_both]
  
  d_obs <- d_state <- d_cond <- NULL
  
  # P(Y)
  if ("observations" %in% type) {
    d_obs <- d[, .(probability = sum(probability)), by = cond_obs]
  }
  # P(Z)
  if ("states" %in% type) {
    d_state <- d[, .(probability = sum(probability)), by = cond_state]
  }
  # P(Y | Z)
  if ("conditionals" %in% type) {
    d_cond <- d[, "probability" := .(probability / sum(probability)), by = cond_state]
  }
  list(observations = d_obs, states = d_state, conditionals = d_cond)
}
#' Predictions from Non-homogeneous Hidden Markov Models
#'
#' @rdname predict
#' @export
predict.nhmm <- function(
    object, newdata, newdata2 = NULL, condition = NULL, 
    type = c("observations", "state", "conditionals"),
    probs = c(0.025, 0.975), boot_idx = TRUE, ...) {
  
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
  }
  if (is.null(condition)) {
    cond <- rep(factor(object$symbol_names), nrow(newdata))
  } else {
    stopifnot_(
      all(condition %in% colnames(newdata)), 
      "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
    )
    cond <- data.frame(
      rep(factor(object$symbol_names, levels = object$symbol_names), nrow(newdata)),
      newdata[rep(seq_len(nrow(newdata)), each = object$n_symbols), condition]
    )
    colnames(cond) <- c(object$channel_names, condition)
  }
  object <- update(object, newdata)
  obs <- create_obsArray(object)
  C <- object$n_channels
  if (C == 1L) {
    out <- predict_nhmm_singlechannel( 
      object$etas$pi, object$X_pi,
      object$etas$A, object$X_A, 
      object$etas$B, object$X_B, 
      array(obs[1, , ], dim(obs)[2:3]),
      object$sequence_lengths, 
      attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
      attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
      attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv")
    )
    if (!is.null(newdata2)) {
      object2 <- update(object, newdata2)
      obs2 <- create_obsArray(object2)
      out2 <- predict_nhmm_singlechannel( 
        object2$etas$pi, object2$X_pi,
        object2$etas$A, object2$X_A, 
        object2$etas$B, object2$X_B, 
        array(obs2[1, , ], dim(obs2)[2:3]),
        object2$sequence_lengths, 
        attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
        attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
        attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv")
      )
    } else {
      out2 <- 0
    }
    mean_estimate <- tapply(stats::na.omit(c(out) - c(out2)), cond, mean)
    
    d <- data.frame(
      expand.grid(attr(mean_estimate, "dimnames")),
      mean = c(mean_estimate)
    )
    if (!is.null(object$boot)) {
      sims <- matrix(NA, nrow(d), length(object$boot_gamma_pi))
      for (i in seq_along(object$boot$gamma_pi)) {
        out <- boot_predict_nhmm_singlechannel( 
          object$etas$pi, object$X_pi,
          object$etas$A, object$X_A, 
          object$etas$B, object$X_B, 
          array(obs[1, , ], dim(obs)[2:3]),
          object$sequence_lengths, 
          attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
          attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
          attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
          object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]]
        )
        if (!is.null(newdata2)) {
          out2 <- boot_predict_nhmm_singlechannel( 
            object2$etas$pi, object2$X_pi,
            object2$etas$A, object2$X_A, 
            object2$etas$B, object2$X_B, 
            array(obs2[1, , ], dim(obs2)[2:3]),
            object2$sequence_lengths, 
            attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
            attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
            attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
            object2$boot$gamma_pi[[i]], object2$boot$gamma_A[[i]], object2$boot$gamma_B[[i]]
          )
        }
        sims[, i] <- c(tapply(na.omit(c(out) - c(out2)), cond, mean))
      }
      
      for(i in seq_along(probs)) {
        d[paste0("q", 100 * probs[i])] <- apply(sims, 1, quantile, probs[i])
      }
    }
  } else {
    stop("Not yet implemented")
  }
  d
}
#' @rdname predict
#' @export
predict.fanhmm <- function(
    object, newdata, newdata2 = NULL, condition = NULL, 
    type = c("observations", "states", "conditionals"),
    probs = c(0.025, 0.975), boot_idx = TRUE, ...) {
  
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
  response_names <- object$channel_names
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
    obs = c(response_names, condition),
    state = c("state", condition),
    both = c("state", response_names, condition)
  )
  
  object <- update(object, newdata)
  obs <- create_obsArray(object, FALSE) # don't set first obs to missing
  stopifnot_(
    object$n_channels == 1L,
    "Multichannel FAN-HMM is not yet supported."
  )
  W <- update_W_for_fanhmm(object, newdata)
  
  out <- unlist(predict_bystate_fanhmm_singlechannel( 
    object$etas$pi, object$X_pi,
    object$etas$A, object$X_A, 
    object$etas$B, object$X_B, 
    array(obs[1, , ], dim(obs)[2:3]),
    object$sequence_lengths, 
    attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
    attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
    attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
    W$W_A, W$W_B, !is.null(object$autoregression_formula)
  ))
  d_mean <- compute_joint(out, newdata, type, conditions, state_names, symbol_names, response_names)
  
  if (!is.null(newdata2)) {
    object2 <- update(object, newdata2)
    obs2 <- create_obsArray(object2, FALSE) # don't set first obs to missing
    W2 <- update_W_for_fanhmm(object2, newdata2)
    out2 <- unlist(predict_bystate_fanhmm_singlechannel( 
      object2$etas$pi, object2$X_pi,
      object2$etas$A, object2$X_A, 
      object2$etas$B, object2$X_B, 
      array(obs2[1, , ], dim(obs2)[2:3]),
      object2$sequence_lengths, 
      attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
      attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
      attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
      W2$W_A, W2$W_B, !is.null(object2$autoregression_formula)
    ))
    d <- compute_joint(
      out2, newdata2, type, conditions, state_names, symbol_names, 
      response_names
    )
    for(i in names(d_mean)) {
      if(!is.null(d_mean[[i]])) {
        d_mean[[i]]$probability <- d_mean[[i]]$probability - d[[i]]$probability
      }
    }
  }
  
  nsim <- length(object$boot$gamma_pi)
  if (nsim > 0) {
    d_boot <- vector("list", nsim)
    boot_idx <- boot_idx & !is.null(object$boot$idx)
    for (i in seq_len(nsim)) {
      out <- simplify2array(boot_predict_bystate_fanhmm_singlechannel( 
        object$etas$pi, object$X_pi,
        object$etas$A, object$X_A, 
        object$etas$B, object$X_B, 
        array(obs[1, , ], dim(obs)[2:3]),
        object$sequence_lengths, 
        attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
        attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
        attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
        W$W_A, W$W_B,
        object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]],
        !is.null(object2$autoregression_formula)
      ))
      if (boot_idx) {
        d_boot[[i]] <- compute_joint(
          out[, , , object$boot$idx[, i], drop = FALSE], 
          newdata, type, conditions, state_names, symbol_names, response_names
        )
      } else {
        d_boot[[i]] <- compute_joint(
          out, newdata, type, conditions, state_names, symbol_names, 
          response_names
        )
      }
      if (!is.null(newdata2)) {
        out2 <- simplify2array(boot_predict_bystate_fanhmm_singlechannel( 
          object2$etas$pi, object2$X_pi,
          object2$etas$A, object2$X_A, 
          object2$etas$B, object2$X_B, 
          array(obs2[1, , ], dim(obs2)[2:3]),
          object2$sequence_lengths, 
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
            response_names
          )
        } else {
          d <- compute_joint(
            out2, newdata2, type, conditions, state_names, symbol_names,
            response_names
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
        qs <- t(apply(di, 2, quantile, probs = probs))
        colnames(qs) <- paste0("q", 100 * probs)
        d_mean[[i]] <- cbind(d_mean[[i]], qs)
      }
    }
  }
  lapply(d_mean[lengths(d_mean) > 0], as.data.frame)
}
#' #' @rdname predict
#' #' @export
#' predict.fanhmm <- function(
    #'     object, newdata, newdata2 = NULL, condition = NULL, 
#'     type = c("observations", "states", "conditionals"),
#'     probs = c(0.025, 0.975), boot_idx = TRUE, ...) {
#'   
#'   type <- match.arg(type, c("observations", "state", "conditionals"), TRUE)
#'   
#'   time <- object$time_variable
#'   id <- object$id_variable
#'   stopifnot_(
#'     !missing(newdata) & is.data.frame(newdata), 
#'     "Argument {.arg newdata} must be a {.cls data.frame} object."
#'   )
#'   stopifnot_(
#'     !is.null(newdata[[id]]), 
#'     "Can't find grouping variable {.var {id}} in {.arg newdata}."
#'   )
#'   stopifnot_(
#'     !is.null(newdata[[time]]), 
#'     "Can't find time index variable {.var {time}} in {.arg newdata}."
#'   )
#'   if (!is.null(newdata2)) {
#'     stopifnot_(
#'       is.data.frame(newdata2), 
#'       "Argument {.arg newdata2} must be a {.cls data.frame} object."
#'     )
#'     stopifnot_(
#'       identical(dim(newdata), dim(newdata2)), 
#'       "Data frames {.arg newdata} and {.arg newdata2} must have the same dimensions."
#'     )
#'     stopifnot_(
#'       identical(names(newdata), names(newdata2)), 
#'       "Data frames {.arg newdata} and {.arg newdata2} must have the same column names."
#'     )
#'     stopifnot_(
#'       identical(newdata[[id]], newdata2[[id]]),
#'       "Grouping variable {.var {id}} must be the same in both data frames."
#'     )
#'     stopifnot_(
#'       identical(newdata[[time]], newdata2[[time]]),
#'       "Time index variable {.var {time}} must be the same in both data frames."
#'     )
#'   }
#'   if (is.null(condition)) {
#'     cond_obs <- data.frame(
#'       y = rep(factor(object$symbol_names), nrow(newdata))
#'     )
#'     names(cond_obs) <- object$channel_names
#'     cond_state <- data.frame(
#'       state = rep(factor(object$state_names), nrow(newdata))
#'     )
#'     cond_cond <- cbind(
#'       state = seq_len(object$n_states), 
#'       cond_obs[rep(seq_len(nrow(cond_obs)), each = object$n_states), ]
#'     )
#'   } else {
#'     stopifnot_(
#'       all(condition %in% colnames(newdata)), 
#'       "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
#'     )
#'     cond_obs <- data.frame(
#'       rep(factor(object$symbol_names, levels = object$symbol_names), nrow(newdata)),
#'       newdata[rep(seq_len(nrow(newdata)), each = object$n_symbols), condition]
#'     )
#'     colnames(cond_obs) <- c(object$channel_names, condition)
#'     cond_state <- data.frame(
#'       rep(factor(object$state_names, levels = object$state_names), nrow(newdata)),
#'       newdata[rep(seq_len(nrow(newdata)), each = object$n_n_states), condition]
#'     )
#'     colnames(cond_state) <- c("state", condition)
#'     cond_cond <- cbind(
#'       state = seq_len(object$n_states), 
#'       cond_obs[rep(seq_len(nrow(cond_obs)), each = object$n_states), ]
#'     )
#'   }
#'   object <- update(object, newdata)
#'   obs <- create_obsArray(object, FALSE) # don't set first obs to missing
#'   stopifnot_(
#'     object$n_channels == 1L,
#'     "Multichannel FAN-HMM is not yet supported."
#'   )
#'   
#'   if (identical(type, "observations")) {
#'     f <- predict_fanhmm_singlechannel
#'     f_boot <- boot_predict_fanhmm_singlechannel
#'   } else {
#'     f <- predict_bystate_fanhmm_singlechannel
#'     f_boot <- boot_predict_bystate_fanhmm_singlechannel
#'   }
#'   
#'   W_A <- W_B <- vector("list", object$n_symbols)
#'   for (i in seq_len(object$n_symbols)) {
#'     d <- newdata
#'     d[[object$channel_names]] <- factor(object$symbol_names[i], levels = object$symbol_names)
#'     mod <- update(object, d)
#'     W_A[[i]] <- mod$X_A
#'     W_B[[i]] <- mod$X_B
#'   }
#'   
#'   out <- unlist(f( 
#'     object$etas$pi, object$X_pi,
#'     object$etas$A, object$X_A, 
#'     object$etas$B, object$X_B, 
#'     array(obs[1, , ], dim(obs)[2:3]),
#'     object$sequence_lengths, 
#'     attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
#'     attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
#'     attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
#'     W_A, W_B, !is.null(object$autoregression_formula)
#'   ))
#'   if (!is.null(newdata2)) {
#'     object2 <- update(object, newdata2)
#'     obs2 <- create_obsArray(object2, FALSE) # don't set first obs to missing
#'     W_A2 <- W_B2 <- vector("list", object2$n_symbols)
#'     for (i in seq_len(object2$n_symbols)) {
#'       d <- newdata2
#'       d[[object2$channel_names]] <- factor(object2$symbol_names[i], levels = object2$symbol_names)
#'       mod <- update(object2, d)
#'       W_A2[[i]] <- mod$X_A
#'       W_B2[[i]] <- mod$X_B
#'     }
#'     out2 <- unlist(f( 
#'       object2$etas$pi, object2$X_pi,
#'       object2$etas$A, object2$X_A, 
#'       object2$etas$B, object2$X_B, 
#'       array(obs2[1, , ], dim(obs2)[2:3]),
#'       object2$sequence_lengths, 
#'       attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
#'       attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
#'       attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
#'       W_A2, W_B2, !is.null(object2$autoregression_formula)
#'     ))
#'   } else {
#'     out2 <- 0
#'   }
#'   
#'   if ("observations" %in% type) {
#'     prob_obs <- tapply(stats::na.omit(c(out - out2)), cond_obs, mean)
#'     d_obs <- data.frame(
#'       expand.grid(attr(prob_obs, "dimnames")),
#'       mean = c(prob_obs)
#'     )
#'   } else {
#'     d_obs <- NULL
#'   }
#'   if ("states" %in% type) {
#'     prob_state <- tapply(stats::na.omit(c(out - out2)), cond_state, mean)
#'     d_state <- data.frame(
#'       expand.grid(attr(prob_state, "dimnames")),
#'       mean = c(prob_state)
#'     )
#'   }
#'   if ("conditionals" %in% type) {
#'     prob_cond <- tapply(stats::na.omit(c(out - out2)), cond_cond, mean)
#'     d_cond <- data.frame(
#'       expand.grid(attr(prob_cond, "dimnames")),
#'       mean = c(prob_cond)
#'     )
#'   }
#'   
#'   nsim <- length(object$boot$gamma_pi)
#'   if (nsim > 0) {
#'     if ("observations" %in% type) {
#'       sims_obs <- matrix(NA, nrow(d), nsim)
#'     }
#'     if ("states" %in% type) {
#'       sims_state <- matrix(NA, nrow(d), nsim)
#'     }
#'     if ("conditionals" %in% type) {
#'       sims_cond <- matrix(NA, nrow(d), nsim)
#'     }
#'     
#'     for (i in seq_len(nsim)) {
#'       out <- f_boot( 
#'         object$etas$pi, object$X_pi,
#'         object$etas$A, object$X_A, 
#'         object$etas$B, object$X_B, 
#'         array(obs[1, , ], dim(obs)[2:3]),
#'         object$sequence_lengths, 
#'         attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
#'         attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
#'         attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
#'         W_A, W_B,
#'         object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]],
#'         !is.null(object2$autoregression_formula)
#'       )
#'       if (!is.null(newdata2)) {
#'         out2 <- f_boot( 
#'           object2$etas$pi, object2$X_pi,
#'           object2$etas$A, object2$X_A, 
#'           object2$etas$B, object2$X_B, 
#'           array(obs2[1, , ], dim(obs2)[2:3]),
#'           object2$sequence_lengths, 
#'           attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
#'           attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
#'           attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
#'           W_A2, W_B2,
#'           object2$boot$gamma_pi[[i]], object2$boot$gamma_A[[i]], object2$boot$gamma_B[[i]],
#'           !is.null(object2$autoregression_formula)
#'         )
#'       }
#'       if (boot_idx & !is.null(object$boot$idx)) {
#'         x <- na.omit(c((out - out2)[, , object$boot$idx[, i], drop = FALSE]))
#'       } else {
#'         x <- na.omit(c(out - out2))
#'       }
#'       if ("observations" %in% type) {
#'         sims_obs[, i] <- c(tapply(x, cond_obs, mean))
#'       }
#'       if ("states" %in% type) {
#'         sims_state[, i] <- c(tapply(x, cond_state, mean))
#'       }
#'       if ("conditionals" %in% type) {
#'         sims_cond[, i] <- c(tapply(x, cond_cond, mean))
#'       }
#'     }
#'     if ("observations" %in% type) {
#'       for(i in seq_along(probs)) {
#'         d_obs[paste0("q", 100 * probs[i])] <- apply(sims_obs, 1, quantile, probs[i])
#'       }
#'     }
#'     if ("states" %in% type) {
#'       for(i in seq_along(probs)) {
#'         d_state[paste0("q", 100 * probs[i])] <- apply(sims_state, 1, quantile, probs[i])
#'       }
#'     }
#'     if ("conditionals" %in% type) {
#'       for(i in seq_along(probs)) {
#'         d_cond[paste0("q", 100 * probs[i])] <- apply(sims_cond, 1, quantile, probs[i])
#'       }
#'     }
#'   }
#'   d
#' }