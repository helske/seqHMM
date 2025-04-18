

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
  ids <- unique(model$data[[model$id_variable]])
  if (io(model$X_pi)) {
    X <- model$X_pi[, 1L, drop = FALSE]
  } else {
    X <- model$X_pi
  }
  S <- model$n_states
  N <- model$n_sequences
  d <- data.table(
    id = rep(ids, each = S),
    state = model$state_names,
    probability = rep_len(c(get_pi_all(model$gammas$gamma_pi, X)), S * N),
    key = "id"
  )
  setnames(d, "id", model$id_variable)
  d
}
#' @rdname initial_probs
#' @export
get_initial_probs.mnhmm <- function(model,...) {
  x <- lapply(split_mnhmm(model), get_initial_probs)
  do.call(rbind, lapply(seq_along(x), \(i) {
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
#' @rdname transition_probs
#' @export
get_transition_probs.nhmm <- function(model, ...) {
  # avoid CRAN check warnings due to NSE
  time <- NULL
  S <- model$n_states
  id_var <- model$id_variable
  time_var <- model$time_variable
  ids <- unique(model$data[[id_var]])
  times <- unique(model$data[[time_var]])
  times <- unlist(lapply(model$sequence_lengths, \(i) times[seq_len(i)]))
  if (!io(model$X_A) && !tv(model$X_A)) {
    X <- model$X_A[, , 1L, drop = FALSE]
  } else {
    X <- model$X_A
  }
  A <- get_A_all(
    model$gammas$gamma_A, X, tv(model$X_A), model$sequence_lengths
  )
  d <- data.table(
    id = rep(ids, S^2 * model$sequence_lengths),
    time = rep(times, each = S^2),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    probability = rep_len(unlist(A), sum(S^2 * model$sequence_lengths)),
    key = c("id", "time")
  )
  d <- d[time != time[1]] # remove A_1
  setnames(d, c("id", "time"), c(id_var, time_var))
  setorderv(d, c(id_var, time_var, "state_from"))
  d
}
#' @rdname transition_probs
#' @export
get_transition_probs.mnhmm <- function(model,...) {
  x <- lapply(split_mnhmm(model), get_transition_probs)
  do.call(rbind, lapply(seq_along(x), \(i) {
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
#' @rdname emission_probs
#' @export
get_emission_probs.nhmm <- function(model, ...) {
  # avoid CRAN check warnings due to NSE
  time <- NULL
  S <- model$n_states
  T_ <- model$length_of_sequences
  C <- model$n_channels
  M <- model$n_symbols
  symbol_names <- model$symbol_names
  id_var <- model$id_variable
  time_var <- model$time_variable
  ids <- unique(model$data[[id_var]])
  times <- unique(model$data[[time_var]])
  times <- unlist(lapply(model$sequence_lengths, \(i) times[seq_len(i)]))

  out <- vector("list", C)
  for (i in seq_len(C)) {
    if (!io(model$X_B[[i]]) && !tv(model$X_B[[i]])) {
      X <- model$X_B[[i]][, , 1L, drop = FALSE]
    } else {
      X <- model$X_B[[i]]
    }
    B <- get_B_all(
      model$gammas$gamma_B[[i]], X, tv(model$X_B[[i]]), model$sequence_lengths
    )
    out[[i]] <- data.table(
      id = rep(ids, S * M[i] * model$sequence_lengths),
      time = rep(times, each = S * M[i]),
      state = model$state_names,
      response = rep(symbol_names[[i]], each = S),
      probability = rep_len(unlist(B), sum(S * M[i] * model$sequence_lengths)),
      key = c("id", "time")
    )
    # if (!is.null(model$autoregression_formula)) {
    #   out[[i]] <- out[[i]][time != time[1]] # remove B_1 as y_1 fixed
    # }
    y <- model$responses[i]
    setnames(out[[i]], c("id", "time", "response"), c(id_var, time_var, y))
    setorderv(out[[i]], c(id_var, time_var, "state"))
  }
  stats::setNames(out, model$responses)
}
#' @rdname emission_probs
#' @export
get_emission_probs.mnhmm <- function(model, ...) {
  x <- lapply(split_mnhmm(model), get_emission_probs)
  do.call(rbind, lapply(seq_along(x), \(i) {
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
get_cluster_probs.mnhmm <- function(model, ...) {
  ids <- unique(model$data[[model$id_variable]])
  if (io(model$X_omega)) {
    X <- model$X_omega[, 1L, drop = FALSE]
  } else {
    X <- model$X_omega
  }
  D <- model$n_clusters
  N <- model$n_sequences
  d <- data.table(
    id = rep(ids, each = D),
    cluster = model$cluster_names,
    probability = rep_len(c(get_omega_all(model$gammas$gamma_omega, X)), D * N),
    key = "id"
  )
  setnames(d, "id", model$id_variable)
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
  data.table(
    cluster = model$cluster_names,
    id = rep(factor(ids), each = model$n_clusters),
    probability = c(t(prior_cluster_probabilities)),
    key = "id"
  )
}
