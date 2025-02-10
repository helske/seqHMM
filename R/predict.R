#' Predictions from Non-homogeneous Hidden Markov Models
#'
#' @rdname predict
#' @export
predict.nhmm <- function(
    object, newdata, newdata2 = NULL, condition = NULL, 
    probs = c(0.025, 0.975), ...) {
  
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
    probs = c(0.025, 0.975), ...) {
  
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
  stopifnot_(
    object$n_channels == 1L,
    "Multichannel FAN-HMM is not yet supported."
  )
  W_A <- W_B <- vector("list", object$n_symbols)
  for (i in seq_len(object$n_symbols)) {
    d <- newdata
    d[[object$channel_names]] <- factor(object$symbol_names[i], levels = object$symbol_names)
    d[[paste0("lag_", object$channel_names)]] <- group_lag(d, id, object$channel_names)
    mod <- update(object, d)
    W_A[[i]] <- mod$X_A
    W_B[[i]] <- mod$X_B
  }
  C <- object$n_channels
  if (C == 1L) {
    out <- predict_fanhmm_singlechannel( 
      object$etas$pi, object$X_pi,
      object$etas$A, object$X_A, 
      object$etas$B, object$X_B, 
      array(obs[1, , ], dim(obs)[2:3]),
      object$sequence_lengths, 
      attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
      attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
      attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
      W_A, W_B
    )
    if (!is.null(newdata2)) {
      object2 <- update(object, newdata2)
      obs2 <- create_obsArray(object2)
      W_A2 <- W_B2 <- vector("list", object2$n_symbols)
      for (i in seq_len(object$n_symbols)) {
        d <- newdata2
        d[[object2$channel_names]] <- factor(object2$symbol_names[i], levels = object2$symbol_names)
        d[[paste0("lag_", object2$channel_names)]] <- group_lag(d, id, object2$channel_names)
        mod <- update(object2, d)
        W_A2[[i]] <- mod$X_A
        W_B2[[i]] <- mod$X_B
      }
      out2 <- predict_fanhmm_singlechannel( 
        object2$etas$pi, object2$X_pi,
        object2$etas$A, object2$X_A, 
        object2$etas$B, object2$X_B, 
        array(obs2[1, , ], dim(obs2)[2:3]),
        object2$sequence_lengths, 
        attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
        attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
        attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
        W_A2, W_B2
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
      sims <- matrix(NA, nrow(d), length(object$boot$gamma_pi))
      for (i in seq_along(object$boot$gamma_pi)) {
        out <- boot_predict_fanhmm_singlechannel( 
          object$etas$pi, object$X_pi,
          object$etas$A, object$X_A, 
          object$etas$B, object$X_B, 
          array(obs[1, , ], dim(obs)[2:3]),
          object$sequence_lengths, 
          attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
          attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
          attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
          W_A, W_B,
          object$boot$gamma_pi[[i]], object$boot$gamma_A[[i]], object$boot$gamma_B[[i]]
        )
        if (!is.null(newdata2)) {
          out2 <- boot_predict_fanhmm_singlechannel( 
            object2$etas$pi, object2$X_pi,
            object2$etas$A, object2$X_A, 
            object2$etas$B, object2$X_B, 
            array(obs2[1, , ], dim(obs2)[2:3]),
            object2$sequence_lengths, 
            attr(object2$X_pi, "icpt_only"), attr(object2$X_A, "icpt_only"), 
            attr(object2$X_B, "icpt_only"), attr(object2$X_A, "iv"), 
            attr(object2$X_B, "iv"), attr(object2$X_A, "tv"), attr(object2$X_B, "tv"),
            W_A2, W_B2,
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