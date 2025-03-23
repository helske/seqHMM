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
#' `NULL`, in which case no quantiles are computed. The quantiles are based on 
#' bootstrap samples of coefficients, stored in `object$boot`.
#'@param condition An optional vector of variable names used for conditional 
#' marginal probabilities. Default is `NULL`, in which case no marginalization is 
#' done.
#' @param newdata An optional data frame containing the new data to be used in 
#' computing the probabilities.
#' @param ... Ignored.
#' @rdname initial_probs
#' @export
get_initial_probs.nhmm <- function(model, probs = NULL, condition = NULL, 
                                   newdata = NULL, ...) {
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  if (!is.null(condition)) {
    if(is.null(newdata)) {
      stopifnot_(
        all(condition %in% colnames(model$data)), 
        "Not all variables defined in {.arg condition} are present in {.arg model$data} ."
      )
    } else {
      stopifnot_(
        all(condition %in% colnames(newdata)), 
        "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
      )
    }
  } else {
    conditions <- "state"
  }
  conditions <- c("state", condition)
  ids <- unique(model$data[[model$id_variable]])
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
    setnames(d, "id", model$id_variable)
  }
  d
}
#' @rdname initial_probs
#' @export
get_initial_probs.mnhmm <- function(model, probs = NULL, condition = NULL, 
                                    newdata = NULL, ...) {
  x <- lapply(split_mnhmm(model), get_initial_probs, probs = probs,
              condition = condition, newdata = newdata)
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
#' @inheritParams get_initial_probs.nhmm
#' @param remove_voids Should the time points corresponding to `TraMineR`'s 
#' void in the observed sequences be removed? Default is `TRUE`.
#' @rdname transition_probs
#' @export
get_transition_probs.nhmm <- function(model, probs = NULL, condition = NULL, 
                                      newdata = NULL, remove_voids = TRUE, ...) {
  # avoid CRAN check warnings due to NSE
  time <- NULL
  stopifnot_(
    checkmate::test_flag(x = remove_voids), 
    "Argument {.arg remove_voids} must be a single {.cls logical} value."
  )
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  if (!is.null(condition)) {
    if(is.null(newdata)) {
      stopifnot_(
        all(condition %in% colnames(model$data)), 
        "Not all variables defined in {.arg condition} are present in {.arg model$data} ."
      )
    } else {
      stopifnot_(
        all(condition %in% colnames(newdata)), 
        "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
      )
    }
  }
  conditions <- c("state", condition)
  S <- model$n_states
  T_ <- model$length_of_sequences
  id_var <- model$id_variable
  time_var <- model$time_variable
  ids <- unique(model$data[[id_var]])
  times <- unique(model$data[[time_var]])
  
  if (!attr(model$X_A, "iv") && !isTRUE(attr(model$W_A, "iv"))) {
    X <- model$X_A[, , 1L, drop = FALSE]
    ids <- "all"
  } else {
    X <- model$X_A
  }
  d <- data.table(
    id = rep(ids, each = S^2 * T_),
    time = rep(times, each = S^2),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = unlist(get_A_all(
      model$gammas$A, X, attr(model$X_A, "tv")
    ))
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
      X, attr(model$X_A, "tv"), probs, model$sequence_lengths
    )
    for(i in seq_along(probs)) {
      d[paste0("q", 100 * probs[i])] <- qs[, i]
    }
  }
  
  d <- d[time == time[1]] # remove A_1
  setnames(d, c("id", "time"), c(id_var, time_var))
  if (remove_voids) {
    stats::na.omit(d)
  } else {
    d
  }
}
#' @rdname transition_probs
#' @export
get_transition_probs.mnhmm <- function(model, probs = NULL, condition = NULL, newdata = NULL, remove_voids = TRUE, ...) {
  x <- lapply(split_mnhmm(model), 
              get_transition_probs, probs = probs, condition = condition, 
              newdata = newdata, remove_voids = remove_voids)
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
#' @inheritParams get_initial_probs.nhmm
#' @param remove_voids Should the time points corresponding to `TraMineR`'s 
#' void in the observed sequences be removed? Default is `TRUE`.
#' @rdname emission_probs
#' @export
get_emission_probs.nhmm <- function(model, probs = NULL, condition = NULL, 
                                    newdata = NULL, remove_voids = TRUE, ...) {
  # avoid CRAN check warnings due to NSE
  time <- NULL
  stopifnot_(
    checkmate::test_flag(x = remove_voids), 
    "Argument {.arg remove_voids} must be a single {.cls logical} value."
  )
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  if (!is.null(condition)) {
    if(is.null(newdata)) {
      stopifnot_(
        all(condition %in% colnames(model$data)), 
        "Not all variables defined in {.arg condition} are present in {.arg model$data} ."
      )
    } else {
      stopifnot_(
        all(condition %in% colnames(newdata)), 
        "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
      )
    }
  }
  conditions <- c("state", condition)
  S <- model$n_states
  C <- model$n_channels
  T_ <- model$length_of_sequences
  M <- model$n_symbols
  ids <- unique(model$data[[model$id_variable]])
  times <- unique(model$data[[model$time_variable]])
  if (C == 1L) {
    symbol_names <- list(model$symbol_names)
    model$gammas$B <- list(model$gammas$B)
  } else {
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
      data.table(
        id = rep(ids, each = S * M[i] * T_),
        time = rep(times, each = S * M[i]),
        state = model$state_names,
        response = model$responses[i],
        observation = rep(symbol_names[[i]], each = S),
        estimate = unlist(get_B_all(
          model$gammas$B[[i]], X, attr(model$X_B, "tv")
        ))
      )
    })
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
        X, attr(model$X_B, "tv"), probs, model$sequence_lengths
      )
    } else {
      qs <- do.call(
        rbind,
        lapply(seq_len(C), function(i) {
          get_B_qs(
            lapply(model$boot$gamma_B, "[[", i), 
            X, attr(model$X_B, "tv"), probs, model$sequence_lengths
          )
        })
      )
    }
    for(i in seq_along(probs)) {
      d[paste0("q", 100 * probs[i])] <- qs[, i]
    }
  }
  if (!is.null(model$autoregression_formula)) {
    # remove B_1 as y_1 is fixed
    d <- d[time == time[1]]
  }
  setnames(d, c("id", "time"), c(model$id_variable, model$time_variable))
  if (remove_voids) {
    stats::na.omit(d)
  } else {
    d
  }
}
#' @rdname emission_probs
#' @export
get_emission_probs.mnhmm <- function(model, probs = NULL, condition = NULL, 
                                     newdata = NULL, remove_voids = TRUE, ...) {
  x <- lapply(split_mnhmm(model), 
              get_emission_probs, probs = probs, condition = condition,
              newdata = newdata, remove_voids = remove_voids)
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
#' @inheritParams get_initial_probs.nhmm
#' @rdname cluster_probs
#' @export
#' @seealso [posterior_cluster_probabilities()].
get_cluster_probs.mnhmm <- function(model, probs = NULL, condition = NULL, 
                                    newdata = NULL, ...) {
  if (!is.null(newdata)) {
    model <- update(model, newdata)
  }
  if (!is.null(condition)) {
    if(is.null(newdata)) {
      stopifnot_(
        all(condition %in% colnames(model$data)), 
        "Not all variables defined in {.arg condition} are present in {.arg model$data} ."
      )
    } else {
      stopifnot_(
        all(condition %in% colnames(newdata)), 
        "Not all variables defined in {.arg condition} are present in {.arg newdata} ."
      )
    }
  }
  conditions <- c("cluster", condition)
  ids <- unique(model$data[[model$id_variable]])
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
  setnames(d, "id", model$id_variable)
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
