#' @rdname initial_probs
#' @export
get_initial_probs <- function(model, ...) {
  UseMethod("get_initial_probs", model)
}
#' @rdname transition_probs
#' @export
get_transition_probs <- function(model, ...) {
  UseMethod("get_transition_probs", model)
}
#' @rdname emission_probs
#' @export
get_emission_probs <- function(model, ...) {
  UseMethod("get_emission_probs", model)
}
#' @rdname cluster_probs
#' @export
get_cluster_probs <- function(model, ...) {
  UseMethod("get_cluster_probs", model)
}
#' Extract the Initial State Probabilities of Hidden Markov Model
#' @param model A hidden Markov model.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.975)`. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param ... Ignored.
#' @rdname initial_probs
#' @export
get_initial_probs.nhmm <- function(model, probs = NULL, 
                                   newdata = NULL, ...) {
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
  }
  if (attr(model$X_pi, "icpt_only")) {
    X <- model$X_pi[, 1L, drop = FALSE]
    ids <- "all"
  } else {
    X <- model$X_pi
  }
  d <- data.frame(
    id = rep(ids, each = model$n_states),
    state = model$state_names,
    estimate = c(
      get_pi_all(model$gammas$pi, X)
    )
  )
  d <- stats::setNames(d, c(model$id_variable, "state", "estimate"))
  if (!is.null(probs)) {
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    B <- length(model$boot$gamma_pi)
    S <- model$n_states
    K <- nrow(X)
    qs <- get_pi_qs(model$boot$gamma_pi, X, probs)
    for(i in seq_along(probs)) {
      d[paste0("q", 100 * probs[i])] <- qs[, i]
    } 
  }
  d
}
#' @rdname initial_probs
#' @export
get_initial_probs.mnhmm <- function(model, probs = NULL, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  x <- lapply(split_mnhmm(model), get_initial_probs, probs = probs)
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(cluster = names(x)[i], x[[i]])
  }))
}
#' @rdname initial_probs
#' @export
get_initial_probs.hmm <- function(model, ...) {
  model$initial_probs
}
#' @rdname initial_probs
#' @export
get_initial_probs.mhmm <- function(model, ...) {
  model$initial_probs
}
#' Extract the State Transition Probabilities of Hidden Markov Model
#' @param model A hidden Markov model.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param remove_voids Should the time points corresponding to `TraMineR`'s 
#' void in the observed sequences be removed? Default is `TRUE`.
#' @param ... Ignored.
#' @rdname transition_probs
#' @export
get_transition_probs.nhmm <- function(model, probs = NULL, newdata = NULL, 
                                      remove_voids = TRUE, ...) {
  stopifnot_(
    checkmate::test_flag(x = remove_voids), 
    "Argument {.arg remove_voids} must be a single {.cls logical} value."
  )
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  S <- model$n_states
  T_ <- model$length_of_sequences
  model$X_A[attr(model$X_A, "missing")] <- NA
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
    times <- as.numeric(colnames(model$observations))
  } else {
    ids <- rownames(model$observations[[1]])
    times <- as.numeric(colnames(model$observations[[1]]))
  }
  if (!attr(model$X_A, "iv") && !isTRUE(attr(model$W_A, "iv"))) {
    X <- model$X_A[, , 1L, drop = FALSE]
    ids <- "all"
  } else {
    X <- model$X_A
  }
  d <- data.frame(
    id = rep(ids, each = S^2 * T_),
    time = rep(times, each = S^2),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = unlist(get_A_all(
      model$gammas$A, X, attr(model$X_A, "tv") || isTRUE(attr(model$W_A, "tv"))
    ))
  )
  d <- stats::setNames(
    d, 
    c(
      model$id_variable, model$time_variable, 
      "state_from", "state_to", "estimate"
    )
  )
  if (!is.null(probs)) {
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    qs <- get_A_qs(
      model$boot$gamma_A, 
      X, attr(model$X_A, "tv") || isTRUE(attr(model$W_A, "tv")), probs
    )
    for(i in seq_along(probs)) {
      d[paste0("q", 100 * probs[i])] <- qs[, i]
    }
  }
  if (remove_voids) {
    d[complete.cases(d), ]
  } else {
    d
  }
}
#' @rdname transition_probs
#' @export
get_transition_probs.mnhmm <- function(model, probs = NULL, newdata = NULL, remove_voids = TRUE, ...) {
  x <- lapply(split_mnhmm(model), 
              get_transition_probs, probs = probs, newdata = newdata, remove_voids = remove_voids)
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(cluster = names(x)[i], x[[i]])
  }))
}
#' @rdname transition_probs
#' @export
get_transition_probs.hmm <- function(model, ...) {
  model$transition_probs
}
#' @rdname transition_probs
#' @export
get_transition_probs.mhmm <- function(model, ...) {
  model$transition_probs
}
#' Extract the Emission Probabilities of Hidden Markov Model
#' @param model A hidden Markov model.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param remove_voids Should the time points corresponding to `TraMineR`'s 
#' void in the observed sequences be removed? Default is `TRUE`.
#' @param ... Ignored.
#' @rdname emission_probs
#' @export
get_emission_probs.nhmm <- function(model, probs = NULL, newdata = NULL, remove_voids = TRUE, ...) {
  stopifnot_(
    checkmate::test_flag(x = remove_voids), 
    "Argument {.arg remove_voids} must be a single {.cls logical} value."
  )
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  S <- model$n_states
  C <- model$n_channels
  T_ <- model$length_of_sequences
  M <- model$n_symbols
  model$X_B[attr(model$X_B, "missing")] <- NA
  if (C == 1L) {
    ids <- rownames(model$observations)
    times <- as.numeric(colnames(model$observations))
    symbol_names <- list(model$symbol_names)
    model$gammas$B <- list(model$gammas$B)
  } else {
    ids <- rownames(model$observations[[1]])
    times <- as.numeric(colnames(model$observations[[1]]))
    symbol_names <- model$symbol_names
  }
  if (!attr(model$X_B, "iv")) {
    X <- model$X_B[, , 1L, drop = FALSE]
    ids <- "all"
  } else {
    X <- model$X_B
  }
  d <- do.call(
    rbind,
    lapply(seq_len(C), function(i) {
      data.frame(
        id = rep(ids, each = S * M[i] * T_),
        time = rep(times, each = S * M[i]),
        state = model$state_names,
        channel = model$channel_names[i],
        observation = rep(symbol_names[[i]], each = S),
        estimate = unlist(get_B_all(
          model$gammas$B[[i]], X, attr(model$X_B, "tv")
        ))
      )
    })
  )
  d <- stats::setNames(
    d, 
    c(model$id_variable, model$time_variable, "state", "channel", 
      "observation", "estimate")
  )
  if (!is.null(probs)) {
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    if (C == 1) {
      qs <- get_B_qs(
        model$boot$gamma_B, 
        X, attr(model$X_B, "tv"), probs
      )
    } else {
      qs <- do.call(
        rbind,
        lapply(seq_len(C), function(i) {
          get_B_qs(
            lapply(model$boot$gamma_B, "[[", i), 
            X, attr(model$X_B, "tv"), probs
          )
        })
      )
    }
    for(i in seq_along(probs)) {
      d[paste0("q", 100 * probs[i])] <- qs[, i]
    }
  }
  if (remove_voids) {
    d[complete.cases(d), ]
  } else {
    d
  }
}
#' @rdname emission_probs
#' @export
get_emission_probs.mnhmm <- function(model, probs = NULL, newdata = NULL, remove_voids = TRUE, ...) {
  x <- lapply(split_mnhmm(model), 
              get_emission_probs, probs = probs, newdata = newdata, remove_voids = remove_voids)
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(cluster = names(x)[i], x[[i]])
  }))
}
#' @rdname emission_probs
#' @export
get_emission_probs.hmm <- function(model, ...) {
  model$emission_probs
}
#' @rdname emission_probs
#' @export
get_emission_probs.mhmm <- function(model, ...) {
  model$emission_probs
}
#' Extract the Prior Cluster Probabilities of MHMM or MNHMM
#' 
#' @param model A mixture hidden Markov model.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param ... Ignored.
#' @rdname cluster_probs
#' @export
#' @seealso [posterior_cluster_probabilities()].
get_cluster_probs.mnhmm <- function(model, probs = NULL, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
  }
  if (attr(model$X_omega, "icpt_only")) {
    X <- model$X_omega[, 1L, drop = FALSE]
    ids <- "all"
  } else {
    X <- model$X_omega
  }
  d <- data.frame(
    cluster = model$cluster_names,
    id = rep(ids, each = model$n_clusters),
    estimate = c(get_omega_all(
      model$gammas$omega, X
    ))
  )
  d <- stats::setNames(d, c("cluster", model$id_variable, "estimate"))
  if (!is.null(probs)) {
    stopifnot_(
      !is.null(model$boot),
      paste0(
        "Model does not contain bootstrap samples of coefficients. ",
        "Run {.fn bootstrap_coefs} first."
      )
    )
    B <- length(model$boot$gamma_omega)
    D <- model$n_clusters
    K <- nrow(X)
    qs <- get_omega_qs(model$boot$gamma_omega, X, probs)
    for(i in seq_along(probs)) {
      d[paste0("q", 100 * probs[i])] <- qs[, i]
    }
  }
  d
}
#' @rdname cluster_probs
#' @export
get_cluster_probs.mhmm <- function(model, ...) {
  pr <- exp(model$X %*% model$coefficients)
  prior_cluster_probabilities <- pr / rowSums(pr)
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
  }
  data.frame(
    cluster = model$cluster_names,
    id = rep(ids, each = model$n_clusters),
    estimate = c(t(prior_cluster_probabilities))
  )
}
#' Get the Estimated Initial, Transition, Emission and (Prior) Cluster 
#' Probabilities for NHMM or MNHMM
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`. The quantiles are based on bootstrap samples of 
#' coefficients, stored in `object$boot`.
#' @param newdata An optional data frame containing the new data to be used in 
#' computing the probabilities.
#' @param remove_voids Should the time points corresponding to `TraMineR`'s 
#' void in the observed sequences be removed? Default is `TRUE`.
#' @param ... Ignored.
#' @rdname get_probs
#' @export
get_probs <- function(model, ...) {
  UseMethod("get_probs", model)
}
#' @rdname get_probs
#' @export
get_probs.nhmm <- function(model, probs = NULL, newdata = NULL, remove_voids = TRUE, ...) {
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  out <- list(
    initial_probs = get_initial_probs(model, probs),
    transition_probs = get_transition_probs(model, probs, remove_voids = remove_voids),
    emission_probs = get_emission_probs(model, probs, remove_voids = remove_voids)
  )
  rownames(out$initial_probs) <- NULL
  rownames(out$transition_probs) <- NULL
  rownames(out$emission_probs) <- NULL
  out
}
#' @rdname get_probs
#' @export
get_probs.mnhmm <- function(model, probs = NULL, newdata = NULL, remove_voids = TRUE, ...) {
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  out <- list(
    initial_probs = get_initial_probs(model, probs),
    transition_probs = get_transition_probs(model, probs, remove_voids = remove_voids),
    emission_probs = get_emission_probs(model, probs, remove_voids = remove_voids),
    cluster_probs = get_cluster_probs(model, probs)
  )
  rownames(out$initial_probs) <- NULL
  rownames(out$transition_probs) <- NULL
  rownames(out$emission_probs) <- NULL
  rownames(out$cluster_probs) <- NULL
  out
}
