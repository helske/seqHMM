#' Predictions from Non-homogeneous Hidden Markov Models
#'
#' @rdname predict
#' @export
predict.nhmm <- function(object, newdata = NULL, ...) {
  
  time <- object$time_variable
  id <- object$id_variable
  if (!is.null(newdata)) {
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
  } else {
    stopifnot_(
      !is.null(object$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    newdata <- object$data
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
    d <- data.frame(
      observation = object$symbol_names,
      time = rep(
        as.numeric(colnames(object$observations)), 
        each = object$n_symbols
      ),
      id = rep(
        rownames(object$observations), 
        each = object$n_symbols * object$length_of_sequences
      ),
      estimate = c(out)
    ) |> 
      stats::na.omit()
    colnames(d)[2] <- time
    if (!is.null(object$boot)) {
      out <- boot_predict_nhmm_singlechannel( 
        object$etas$pi, object$X_pi,
        object$etas$A, object$X_A, 
        object$etas$B, object$X_B, 
        array(obs[1, , ], dim(obs)[2:3]),
        object$sequence_lengths, 
        attr(object$X_pi, "icpt_only"), attr(object$X_A, "icpt_only"), 
        attr(object$X_B, "icpt_only"), attr(object$X_A, "iv"), 
        attr(object$X_B, "iv"), attr(object$X_A, "tv"), attr(object$X_B, "tv"),
        object$boot$gamma_pi, object$boot$gamma_A, object$boot$gamma_B
      )
      d_boot <- data.frame(
        observation = object$symbol_names,
        time = rep(
          as.numeric(colnames(object$observations)), 
          each = object$n_symbols
        ),
        id = rep(
          rownames(object$observations), 
          each = object$n_symbols * object$length_of_sequences
        ),
        estimate = unlist(out$obs_prob)
      ) |> 
        stats::na.omit()
      colnames(d_boot)[2] <- time
      d$type <- "MLE"
      d_boot$type <- "Bootstrap"
      d <- rbind(d, d_boot)
    }
  } else {
    stop("Not yet implemented")
  }
  d
}

#' @rdname predict
#' @export
predict.fanhmm <- function(object, newdata = NULL, ...) {
  
  time <- object$time_variable
  id <- object$id_variable
  if (!is.null(newdata)) {
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
  } else {
    stopifnot_(
      !is.null(object$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    newdata <- object$data
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
    W_A[[i]] <- mod$X_A[, , , drop = FALSE]
    W_B[[i]] <- mod$X_B[, , , drop = FALSE]
  }
  
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
  d <- data.frame(
    observation = object$symbol_names,
    time = rep(
      as.numeric(colnames(object$observations)), 
      each = object$n_symbols
    ),  
    id = rep(
      rownames(object$observations), 
      each = object$n_symbols * object$length_of_sequences
    ),
    estimate = c(out)
  ) |> 
    stats::na.omit()
  colnames(d)[2] <- time
  
  if (!is.null(object$boot)) {
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
      object$boot$gamma_pi, object$boot$gamma_A, object$boot$gamma_B
    )
    d_boot <- data.frame(
      observation = object$symbol_names,
      time = rep(
        as.numeric(colnames(object$observations)), 
        each = object$n_symbols
      ),
      id = rep(
        rownames(object$observations), 
        each = object$n_symbols * object$length_of_sequences
      ),
      estimate = unlist(out$obs_prob)
    ) |> 
      stats::na.omit()
    colnames(d_boot)[2] <- time
    d$type <- "MLE"
    d_boot$type <- "Bootstrap"
    d <- rbind(d, d_boot)
  }
  d
}