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
#' @param ... Ignored.
#' @rdname initial_probs
#' @export
get_initial_probs.nhmm <- function(model, ...) {
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
  }
  if (!attr(model, "iv_pi")) {
    d <- data.frame(
      id = rep(ids, each = model$n_states),
      state = model$state_names,
      estimate = get_pi(
        model$coefficients$gamma_pi_raw, model$X_initial[, 1L], FALSE
      )
    )
  } else {
    d <- data.frame(
      id = rep(ids, each = model$n_states),
      state = model$state_names,
      estimate = c(
        get_pi_all(model$coefficients$gamma_pi_raw, model$X_initial, FALSE)
      )
    )
  }
  stats::setNames(d, c(model$id_variable, "state", "estimate"))
}
#' @rdname initial_probs
#' @export
get_initial_probs.mnhmm <- function(model, ...) {
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
  }
  if (!attr(model, "iv_pi")) {
    d <- do.call(
      rbind,
      lapply(seq_len(model$n_clusters), function(i) {
        data.frame(
          cluster = model$cluster_names[i],
          id = rep(ids, each = model$n_states),
          state = model$state_names[[i]],
          estimate = c(
            get_pi(model$coefficients$gamma_pi_raw[[i]], model$X_initial[, 1L], 
                   FALSE)
          )
        )
      })
    )
  } else {
    d <- do.call(
      rbind,
      lapply(seq_len(model$n_clusters), function(i) {
        data.frame(
          cluster = model$cluster_names[i],
          id = rep(ids, each = model$n_states),
          state = model$state_names[[i]],
          estimate = c(
            get_pi_all(
              model$coefficients$gamma_pi_raw[[i]], model$X_initial, FALSE
            )
          )
        )
      })
    )
  }
  stats::setNames(d, c("cluster", model$id_variable, "state", "estimate"))
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
#' @param ... Ignored.
#' @rdname transition_probs
#' @export
get_transition_probs.nhmm <- function(model, ...) {
  S <- model$n_states
  T_ <- model$length_of_sequences
  model$X_transition[attr(model, "missing_X_transition")] <- NA
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
    times <- colnames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
    times <- colnames(model$observations[[1]])
  }
  if (!attr(model, "iv_A")) {
    X <- matrix(model$X_transition[, , 1L], ncol = model$length_of_sequences)
    d <- data.frame(
      id = rep(ids, each = S^2 * T_),
      time = rep(times, each = S^2),
      state_from = model$state_names,
      state_to = rep(model$state_names, each = S),
      estimate = c(get_A(
        model$coefficients$gamma_A_raw, X, FALSE, attr(model, "tv_A")
      ))
    )
  } else {
    d <- data.frame(
      id = rep(ids, each = S^2 * T_),
      time = rep(times, each = S^2),
      state_from = model$state_names,
      state_to = rep(model$state_names, each = S),
      estimate = unlist(
        get_A_all(
          model$coefficients$gamma_A_raw, model$X_transition, FALSE, 
          attr(model, "tv_A")
        )
      )
    )
  }
  stats::setNames(
    d, 
    c(
      model$id_variable, model$time_variable, 
      "state_from", "state_to", "estimate"
    )
  )
}
#' @rdname transition_probs
#' @export
get_transition_probs.mnhmm <- function(model, ...) {
  S <- model$n_states
  T_ <- model$length_of_sequences
  D <- model$n_clusters
  model$X_transition[attr(model, "missing_X_transition")] <- NA
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
    times <- colnames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
    times <- colnames(model$observations[[1]])
  }
  if (!attr(model, "iv_A")) {
    X <- matrix(model$X_transition[, , 1L], ncol = model$length_of_sequences)
    d <- do.call(
      rbind,
      lapply(seq_len(D), function(i) {
        data.frame(
          cluster = model$cluster_names[i],
          id = rep(ids, each = S^2 * T_),
          time = rep(times, each = S^2),
          state_from = model$state_names[[i]],
          state_to = rep(model$state_names[[i]], each = S),
          estimate = c(get_A(
            model$coefficients$gamma_A_raw[[i]], X, FALSE, attr(model, "tv_A")
          ))
        )
      })
    )
  } else {
    d <- do.call(
      rbind,
      lapply(seq_len(D), function(i) {
        data.frame(
          cluster = model$cluster_names[i],
          id = rep(ids, each = S^2 * T_),
          time = rep(times, each = S^2),
          state_from = model$state_names[[i]],
          state_to = rep(model$state_names[[i]], each = S),
          estimate = unlist(
            get_A_all(
              model$coefficients$gamma_A_raw[[i]], model$X_transition, FALSE, 
              attr(model, "tv_A")
            )
          )
        )
      })
    )
  }
  stats::setNames(
    d, 
    c("cluster", model$id_variable, model$time_variable, 
      "state_from", "state_to", "estimate")
  )
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
#' @param ... Ignored.
#' @rdname emission_probs
#' @export
get_emission_probs.nhmm <- function(model, ...) {
  S <- model$n_states
  C <- model$n_channels
  T_ <- model$length_of_sequences
  M <- model$n_symbols
  model$X_emission[attr(model, "missing_X_emission")] <- NA
  if (C == 1L) {
    ids <- rownames(model$observations)
    times <- colnames(model$observations)
    symbol_names <- list(model$symbol_names)
    model$coefficients$gamma_B_raw <- list(model$coefficients$gamma_B_raw)
  } else {
    ids <- rownames(model$observations[[1]])
    times <- colnames(model$observations[[1]])
    symbol_names <- model$symbol_names
  }
  if (!attr(model, "iv_B")) {
    X <- matrix(model$X_emission[, , 1L], ncol = model$length_of_sequences)
    d <- do.call(
      rbind,
      lapply(seq_len(C), function(i) {
        data.frame(
          id = rep(ids, each = S * M[i] * T_),
          time = rep(times, each = S * M[i]),
          state = model$state_names,
          channel = model$channel_names[i],
          observation = rep(symbol_names[[i]], each = S),
          estimate = c(get_B(
            model$coefficients$gamma_B_raw[[i]], X, FALSE, FALSE, 
            attr(model, "tv_B"))
          )
        )
      })
    )
  } else {
    d <- do.call(
      rbind,
      lapply(seq_len(C), function(i) {
        data.frame(
          id = rep(ids, each = S * M[i] * T_),
          time = rep(times, each = S * M[i]),
          state = model$state_names,
          channel = model$channel_names[i],
          observation = rep(symbol_names[[i]], each = S),
          estimate = unlist(
            get_B_all(
              model$coefficients$gamma_B_raw[[i]], model$X_emission, 
              FALSE, FALSE, attr(model, "tv_B")
            )
          )
        )
      })
    )
  }
  stats::setNames(
    d, 
    c(model$id_variable, model$time_variable, "state", "channel", 
      "observation", "estimate")
  )
}
#' @rdname emission_probs
#' @export
get_emission_probs.mnhmm <- function(model, ...) {
  S <- model$n_states
  C <- model$n_channels
  D <- model$n_clusters
  T_ <- model$length_of_sequences
  M <- model$n_symbols
  model$X_emission[attr(model, "missing_X_emission")] <- NA
  if (C == 1L) {
    ids <- rownames(model$observations)
    times <- colnames(model$observations)
    symbol_names <- list(model$symbol_names)
    model$coefficients$gamma_B_raw <- lapply(
      model$coefficients$gamma_B_raw, list
    )
  } else {
    ids <- rownames(model$observations[[1]])
    times <- colnames(model$observations[[1]])
    symbol_names <- model$symbol_names
  }
  if (!attr(model, "iv_B")) {
    X <- matrix(model$X_emission[, , 1L], ncol = model$length_of_sequences)
    d <- do.call(
      rbind,
      lapply(seq_len(D), function(j) {
        do.call(
          rbind,
          lapply(seq_len(C), function(i) {
            data.frame(
              cluster = model$cluster_names[j],
              id = rep(ids, each = S * M[i] * T_),
              time = rep(times, each = S * M[i]),
              state = model$state_names[[j]],
              channel = model$channel_names[i],
              observation = rep(symbol_names[[i]], each = S),
              estimate = c(get_B(
                model$coefficients$gamma_B_raw[[j]][[i]], X, FALSE, FALSE, 
                attr(model, "tv_B"))
              )
            )
          })
        )
      })
    )
  } else {
    d <- do.call(
      rbind,
      lapply(seq_len(D), function(j) {
        do.call(
          rbind,
          lapply(seq_len(C), function(i) {
            data.frame(
              cluster = model$cluster_names[j],
              id = rep(ids, each = S * M[i] * T_),
              time = rep(times, each = S * M[i]),
              state = model$state_names[[j]],
              channel = model$channel_names[i],
              observation = rep(symbol_names[[i]], each = S),
              estimate = unlist(
                get_B_all(
                  model$coefficients$gamma_B_raw[[j]][[i]], model$X_emission, 
                  FALSE, FALSE,attr(model, "tv_B")
                )
              )
            )
          })
        )
      })
    )
  }
  stats::setNames(
    d, 
    c("cluster", model$id_variable, model$time_variable, "state", "channel", 
      "observation", "estimate")
  )
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
#' @param ... Ignored.
#' @rdname cluster_probs
#' @export
#' @seealso [posterior_cluster_probabilities()].
get_cluster_probs.mnhmm <- function(model, ...) {
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
  }
  if (!attr(model, "iv_omega")) {
    d <- data.frame(
      cluster = model$cluster_names,
      id = rep(ids, each = model$n_clusters),
      estimate = c(get_omega(
        model$coefficients$gamma_omega_raw, model$X_cluster[, 1L], FALSE
      ))
    )
  } else {
    d <- data.frame(
      cluster = model$cluster_names,
      id = rep(ids, each = model$n_clusters),
      estimate = c(
        get_omega_all(
          model$coefficients$gamma_omega_raw, model$X_cluster, FALSE
        )
      )
    )
  }
  stats::setNames(d, c("cluster", model$id_variable, "estimate"))
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
get_probs.nhmm <- function(model, newdata = NULL, remove_voids = TRUE, ...) {
  stopifnot_(
    checkmate::test_flag(x = remove_voids), 
    "Argument {.arg remove_voids} must be a single {.cls logical} value."
  )
  if (!is.null(newdata)) {
    time <- model$time_variable
    id <- model$id_variable
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
    object <- update(object, newdata = newdata)
  }
  S <- model$n_states
  M <- model$n_symbols
  C <- model$n_channels
  N <- model$n_sequences
  T_ <- model$length_of_sequences
  out <- list(
    initial_probs = get_initial_probs(model),
    transition_probs = get_transition_probs(model),
    emission_probs = get_emission_probs(model)
  )
  rownames(out$initial_probs) <- NULL
  rownames(out$transition_probs) <- NULL
  rownames(out$emission_probs) <- NULL
  if (remove_voids) {
    list(
      initial_probs = out$initial_probs, 
      transition_probs = out$transition_probs[complete.cases(out$transition_probs), ],
      emission_probs = out$emission_probs[complete.cases(out$emission_probs), ]
    )
  } else out
}
#' @rdname get_probs
#' @export
get_probs.mnhmm <- function(model, newdata = NULL, remove_voids = TRUE, ...) {
  
  stopifnot_(
    checkmate::test_flag(x = remove_voids), 
    "Argument {.arg remove_voids} must be a single {.cls logical} value."
  )
  if (!is.null(newdata)) {
    time <- model$time_variable
    id <- model$id_variable
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
    object <- update(object, newdata = newdata)
  }
  S <- model$n_states
  M <- model$n_symbols
  C <- model$n_channels
  N <- model$n_sequences
  T_ <- model$length_of_sequences
  out <- list(
    initial_probs = get_initial_probs(model),
    transition_probs = get_transition_probs(model),
    emission_probs = get_emission_probs(model),
    cluster_probs = get_cluster_probs(model)
  )
  rownames(out$initial_probs) <- NULL
  rownames(out$transition_probs) <- NULL
  rownames(out$emission_probs) <- NULL
  rownames(out$cluster_probs) <- NULL
  if (remove_voids) {
    list(
      initial_probs = out$initial_probs, 
      transition_probs = out$transition_probs[complete.cases(out$transition_probs), ],
      emission_probs = out$emission_probs[complete.cases(out$emission_probs), ],
      cluster_probs = out$cluster_probs
    )
  } else out
}
