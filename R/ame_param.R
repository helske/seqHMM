#' Average Marginal Effects on NHMM Parameters
#' 
#' The function `ame_param` computes the average marginal effect (AME) of the 
#' model covariate \eqn{X} on the model parameters by marginalizing over the sequences. 
#' Under the assumption of no unobserved confounding (i.e., there are no 
#' unobserved variables that influence the covariate \eqn{X} and the outcome), 
#' these can be regarded as the causal effects of the covariate on 
#' the initial,  emission, and transition probabilities of the model. In case 
#' `values` argument is a single value \eqn{x}, the function returns the 
#' interventional initial, transition, emission probabilities
#' \deqn{P(z_1 | do(X_1 = x))}
#' \deqn{P(z_t | do(X_{t-1} = x), z_{t-1})}
#' \deqn{P(y_t | do(X_t = x), z_t)}
#' and in a case `values` contains two values \eqn{x} and \eqn{w} a shift in 
#' interventional distributions, i.e.,
#' \deqn{P(z_1 | do(X_1 = x)) - P(z_1 | do(X_1 = w))}
#' \deqn{P(z_t | do(X_{t-1} = x), z_{t-1}) - P(z_t | do(X_{t-1} = w), z_{t-1})}
#' \deqn{P(y_t | do(X_t = x), z_t) - P(y_t | do(X_t = w), z_t)}.
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing one or two values for `variable`. 
#' See details.
#' @param newdata Optional data frame which is used for marginalization.
#' @param probs Quantiles of interest of average marginal effect.
#' @param ... Ignored.
#' @rdname ame_param
#' @export
ame_param <- function(model, variable, values, ...) {
  UseMethod("ame_param", model)
}
#' @rdname ame_param
#' @export
ame_param.nhmm <- function(
    model, variable, values, newdata = NULL, probs, ...) {
  stopifnot_(
    attr(model, "intercept_only") == FALSE,
    "Model does not contain any covariates."
  )
  stopifnot_(
    checkmate::test_string(x = variable), 
    "Argument {.arg variable} must be a single character string."
  )
  stopifnot_(
    length(values) == 2, 
    "Argument {.arg values} should contain two values for 
    variable {.var variable}.")
  
  if (!missing(probs)) {
    return_quantiles <- TRUE
    stopifnot_(
      checkmate::test_numeric(
        x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L
      ),
      "Argument {.arg probs} must be a {.cls numeric} vector with values
      between 0 and 1."
    )
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first in order to compute quantiles."
      )
    )
  } else {
    return_quantiles <- FALSE
  }
  
  time <- model$time_variable
  id <- model$id_variable
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
    stopifnot_(
      !is.null(newdata[[variable]]), 
      "Can't find time variable {.var {variable}} in {.arg newdata}."
    )
  } else {
    stopifnot_(
      !is.null(model$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    newdata <- model$data
  }
  newdata[[variable]][] <- values[1]
  model1 <- update(model, newdata)
  newdata[[variable]][] <- values[2]
  model2 <- update(model, newdata)
  C <- model$n_channels
  if (C == 1L) {
    times <- as.numeric(colnames(model$observations))
    symbol_names <- list(model$symbol_names)
    model$gammas$B <- list(model$gammas$B)
  } else {
    times <- as.numeric(colnames(model$observations[[1]]))
    symbol_names <- model$symbol_names
  }
  if (attr(model$X_pi, "icpt_only")) {
    X1 <- model1$X_pi[, 1L, drop = FALSE]
    X2 <- model2$X_pi[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_pi
    X2 <- model2$X_pi
  }
  ame_pi <-  data.frame(
    state = model$state_names,
    estimate = rowMeans(
      get_pi_all(model$gammas$pi, X1) - get_pi_all(model$gammas$pi, X2) 
    )
  )
  if (return_quantiles) {
    qs_pi <- get_pi_ame(model$boot$gamma_pi, X1, X2, probs)
    colnames(qs_pi) <- paste0("q", 100 * probs)
    ame_pi <- cbind(ame_pi, qs_pi)
  }
  
  model1$X_A[attr(model$X_A, "missing")] <- NA
  model2$X_A[attr(model$X_A, "missing")] <- NA
  if (!attr(model$X_A, "iv")) {
    X1 <- model1$X_A[, 1L, , drop = FALSE]
    X2 <- model2$X_A[, 1L, , drop = FALSE]
  } else {
    X1 <- model1$X_A
    X2 <- model2$X_A
  }
  tv_A <- attr(model$X_A, "tv")
  S <- model$n_states
  N <- model$n_sequences
  T_ <- model$length_of_sequences
  ame_A <- data.frame(
    time = rep(times, each = S^2),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = c(
      apply(
        array(
          unlist(get_A_all(model$gammas$A, X1, tv_A)) -
            unlist(get_A_all(model$gammas$A, X2, tv_A)), 
          c(S, S, T_, N)
        ), 
        1:3, mean)
    )
  )
  if (return_quantiles) {
    qs_A <- get_A_ame(model$boot$gamma_A, X1, X2, tv_A, probs)
    colnames(qs_A) <- paste0("q", 100 * probs)
    ame_A <- cbind(ame_A, qs_A)
  }
  colnames(ame_A)[1] <- model$time_variable
  
  model1$X_B[attr(model$X_B, "missing")] <- NA
  model1$X_B[attr(model$X_B, "missing")] <- NA
  if (!attr(model$X_B, "iv")) {
    X1 <- model1$X_B[, 1L, drop = FALSE]
    X2 <- model2$X_B[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_B
    X2 <- model2$X_B
  }
  tv_B <- attr(model$X_B, "tv")
  M <- model$n_symbols
  ame_B <- do.call(
    rbind,
    lapply(seq_len(C), function(i) {
      data.frame(
        time = rep(times, each = S * M[i]),
        state = model$state_names,
        channel = model$channel_names[i],
        observation = rep(symbol_names[[i]], each = S),
        estimate = c(
          apply(
            array(
              unlist(get_B_all(model$gammas$B[[i]], X1, tv_B)) -
                unlist(get_B_all(model$gammas$B[[i]], X2, tv_B)), 
              c(S, M[i], T_, N)
            ), 
            1:3, mean)
        )
      )
    })
  )
  if (return_quantiles) {
    if (C == 1) {
      qs_B <- get_B_ame(
        model$boot$gamma_B, X1, X2, tv_B, probs
      )
    } else {
      qs_B <- do.call(
        rbind,
        lapply(seq_len(C), function(i) {
          get_B_ame(
            lapply(model$boot$gamma_B, "[[", i), 
            X1, X2, tv_B, probs
          )
        })
      )
    }
    colnames(qs_B) <- paste0("q", 100 * probs)
    ame_B <- cbind(ame_B, qs_B)
  }
  colnames(ame_B)[1] <- model$time_variable
  
  out <- list(
    initial = ame_pi,
    transition = ame_A,
    emission = ame_B
  )
  class(out) <- "ame_param"
  attr(out, "model") <- "nhmm"
  out
}

#' @rdname ame_param
#' @export
ame_param.mnhmm <- function(
    model, variable, values, newdata = NULL, probs,
    ...) {
  
  x <- lapply(
    split_mnhmm(model), ame_param, variable = variable, values = values, 
    newdata = newdata, probs = probs
  )
  out <- lapply(c("pi", "A", "B"), function(z) {
    do.call(rbind, lapply(seq_along(x), function(i) {
      cbind(cluster = names(x)[i], x[[i]][[z]])
    }))
  }
  )
  names(out) <- c("initial", "transition", "emission")
  
  # no need for checks, ame.nhmm already checks these
  if (is.null(newdata)) {
    newdata <- model$data
  }
  newdata[[variable]][] <- values[1]
  model1 <- update(model, newdata)
  newdata[[variable]][] <- values[2]
  model2 <- update(model, newdata)
  
  if (!attr(model$X_omega, "icpt_only")) {
    X1 <- model1$X_omega[, 1L, drop = FALSE]
    X2 <- model2$X_omega[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_omega
    X2 <- model2$X_omega
  }
  out$cluster <- data.frame(
    cluster = model$cluster_names,
    estimate = rowMeans(
      get_omega_all(model$gammas$omega, X1) - 
        get_omega_all(model$gammas$omega, X2) 
    )
  )
  if (!missing(probs) && !is.null(model$boot)) {
    qs_omega <- get_omega_ame(model$boot$gamma_omega, X1, X2, probs)
    colnames(qs_omega) <- paste0("q", 100 * probs)
    out$cluster <- cbind(out$cluster, qs_omega)
  }
  class(out) <- "ame_param"
  attr(out, "model") <- "mnhmm"
  out
}
