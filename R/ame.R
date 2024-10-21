#' Average Marginal Effects for Non-homogenous Hidden Markov Models
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing one or two values for `variable`.
#' @param newdata Optional data frame which is used for marginalization.
#' @param probs Quantiles of interest of average marginal effect.
#' @param ... Ignored.
#' @rdname ame
#' @export
ame <- function(model, variable, values, ...) {
  UseMethod("ame", model)
}
#' @rdname ame
#' @export
ame.nhmm <- function(
    model, variable, values, newdata = NULL, probs = c(0.05, 0.95),
    ...) {
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
  stopifnot_(
    !is.null(model$boot),
    paste0(
      "Model does not contain bootstrap samples of coefficients. ",
      "Run {.fn bootstrap_coefs} first."
    )
  )
  newdata[[variable]] <- values[1]
  model1 <- update(model, newdata)
  newdata[[variable]] <- values[2]
  model2 <- update(model, newdata)
  C <- model$n_channels
  if (C == 1L) {
    times <- colnames(model$observations)
    symbol_names <- list(model$symbol_names)
  } else {
    times <- colnames(model$observations[[1]])
    symbol_names <- model$symbol_names
  }
  if (!attr(model$X_initial, "iv")) {
    X1 <- model1$X_initial[, 1L, drop = FALSE]
    X2 <- model2$X_initial[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_initial
    X2 <- model2$X_initial
  }
  qs_pi <- get_pi_ame(model$boot$gamma_pi, X1, X2, probs)
  colnames(qs_pi) <- paste0("q", 100 * probs)
  ame_pi <- cbind(
    data.frame(
      state = model$state_names,
      estimate = rowMeans(
        get_pi_all(model$gammas$pi, X1) - get_pi_all(model$gammas$pi, X2) 
      )
    ),
    qs_pi
  )
  model1$X_transition[attr(model$X_transition, "missing")] <- NA
  model2$X_transition[attr(model$X_transition, "missing")] <- NA
  if (!attr(model$X_transition, "iv")) {
    X1 <- model1$X_transition[, 1L, , drop = FALSE]
    X2 <- model2$X_transition[, 1L, , drop = FALSE]
  } else {
    X1 <- model1$X_transition
    X2 <- model2$X_transition
  }
  tv_A <- attr(model$X_transition, "tv")
  S <- model$n_states
  N <- model$n_sequences
  T_ <- model$length_of_sequences
  qs_A <- get_A_ame(model$boot$gamma_A, X1, X2, tv_A, probs)
  colnames(qs_A) <- paste0("q", 100 * probs)
  ame_A <- cbind(
    data.frame(
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
    ),
    qs_A
  )
  colnames(ame_A)[1] <- model$time_variable
  
  model1$X_emission[attr(model$X_emission, "missing")] <- NA
  model1$X_emission[attr(model$X_emission, "missing")] <- NA
  if (!attr(model$X_emission, "iv")) {
    X1 <- model1$X_emission[, 1L, drop = FALSE]
    X2 <- model2$X_emission[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_emission
    X2 <- model2$X_emission
  }
  tv_B <- attr(model$X_emission, "tv")
  M <- model$n_symbols
  if (C == 1) {
    qs_B <- get_B_ame(
      model$boot$gamma_B, X1, X2, tv_B, probs
    )
    model$gammas$B <- list(model$gammas$B)
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
  ame_B <- cbind(
    do.call(
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
    ),
    qs_B
  )
  colnames(ame_B)[1] <- model$time_variable
  
  out <- list(
    initial = ame_pi,
    transition = ame_A,
    emission = ame_B
  )
  class(out) <- "amp"
  attr(out, "model") <- "nhmm"
  out
}

#' @rdname ame
#' @export
ame.mnhmm <- function(
    model, variable, values, newdata = NULL, probs = c(0.05, 0.95),
    ...) {
  
  x <- lapply(
    split_mnhmm(model), ame, variable = variable, values = values, 
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
  newdata[[variable]] <- values[1]
  model1 <- update(model, newdata)
  newdata[[variable]] <- values[2]
  model2 <- update(model, newdata)
  
  if (!attr(model$X_omega, "iv")) {
    X1 <- model1$X_cluster[, 1L, drop = FALSE]
    X2 <- model2$X_cluster[, 1L, drop = FALSE]
  } else {
    X1 <- model1$X_cluster
    X2 <- model2$X_cluster
  }
  qs_omega <- get_omega_ame(model$boot$gamma_omega, X1, X2, probs)
  colnames(qs_omega) <- paste0("q", 100 * probs)
  out$cluster <- cbind(
    data.frame(
      cluster = model$cluster_names,
      estimate = rowMeans(
        get_omega_all(model$gammas$omega, X1) - 
          get_omega_all(model$gammas$omega, X2) 
      )
    ),
    qs_omega
  )
  class(out) <- "amp"
  attr(out, "model") <- "mnhmm"
  out
}
