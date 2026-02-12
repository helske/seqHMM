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
#' @param ids Optional vector of sequence IDs. If `NULL` (the default), 
#' probabilities are computed for all sequences in the data, unless the 
#' probabilities are equal for all sequences, i.e., covariates do not depend on 
#' ID.
#' @return A `data.table` containing the initial state probabilities for each
#' sequence defined by `ids`.
#' @rdname initial_probs
#' @export
get_initial_probs.nhmm <- function(model, ids = NULL, ...) {
  if (!is.null(ids)) {
    missing_ids <- setdiff(ids, unique(model$data[[id]]))
    stopifnot_(
      length(missing_ids) == 0L,
      "The following IDs are not present in the model data: 
      {paste(missing_ids, collapse = ', ')}."
    )
    ids <- factor(ids, levels = levels(model$data[[model$id_variable]]))
  } else {
    ids <- unique(model$data[[model$id_variable]])
  }
  if (io(model$X_pi)) {
    X <- model$X_pi[, 1L, drop = FALSE]
  } else {
    X <- model$X_pi[, ids]
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
  d[]
}
#' @rdname initial_probs
#' @export
get_initial_probs.mnhmm <- function(model, ids = NULL, ...) {
  x <- lapply(split_mnhmm(model), get_initial_probs, ids = ids)
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
get_transition_probs.nhmm <- function(model, ids = NULL, ...) {
  probability <- state_from <- state_to <- S2 <- NULL
  S <- model$n_states
  id <- model$id_variable
  time <- model$time_variable
  states <- model$state_names
  
  if (!is.null(ids)) {
    missing_ids <- setdiff(ids, unique(model$data[[id]]))
    stopifnot_(
      length(missing_ids) == 0L,
      "The following IDs are not present in the model data: 
      {paste(missing_ids, collapse = ', ')}."
    )
    ids <- factor(ids, levels = levels(model$data[[id]]))
  } else {
    ids <- unique(model$data[[id]])
  }
  d <- model$data[list(ids), nomatch = 0L]
  d <- d[rep(seq_row(d), each = S2), 
         list(id, time), 
         env = list(id = id, time = time, S2 = S^2)]
  n <- nrow(d)
  set(d, j = "state_from", value = rep_len(states, n))
  set(d, j = "state_to", value = rep_len(rep(states, each = S), n))
  
  if (!iv(model$X_A)) {
    A <- get_A_all(
      model$gammas$gamma_A, model$X_A[1], tv(model$X_A)
    )
    set(d, j = "probability", value = rep_len(unlist(A), n))
  } else {
    A <- get_A_all(model$gammas$gamma_A, model$X_A[ids], tv(model$X_A))
    set(d, j = "probability", value = unlist(A))
  }
  setorderv(d, c(id, time, "state_from"))
  d[]
}
#' @rdname transition_probs
#' @export
get_transition_probs.mnhmm <- function(model, ids = NULL, ...) {
  x <- lapply(
    split_mnhmm(model), 
    get_transition_probs, ids = ids
  )
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
#' @inheritParams get_transition_probs.nhmm
#' @rdname emission_probs
#' @export
get_emission_probs.nhmm <- function(model, ids = NULL, ...) {
  probability <- state <- SM <- NULL
  responses <- model$responses
  M <- model$n_symbols
  S <- model$n_states
  id <- model$id_variable
  time <- model$time_variable
  states <- model$state_names
  symbols <- model$symbol_names
  out <- stats::setNames(vector("list", length(responses)), responses)
  
  if (!is.null(ids)) {
    missing_ids <- setdiff(ids, unique(model$data[[id]]))
    stopifnot_(
      length(missing_ids) == 0L,
      "The following IDs are not present in the model data: 
      {paste(missing_ids, collapse = ', ')}."
    )
    ids <- factor(ids, levels = levels(model$data[[id]]))
  } else {
    ids <- unique(model$data[[id]])
  }
  if (inherits(model, "fanhmm") && !identical(model$prior_obs, 0L)) {
    B1 <- get_B1(
      model$gammas$gamma_B, model$n_symbols, model$W_X_B, model$prior_obs
    )
    id_idx <- match(ids, unique(model$data[[id]]))
  }
  for (i in seq_along(responses)) {
    y <- responses[i]
    d <- model$data[list(ids), nomatch = 0L]
    d <- d[rep(seq_row(d), each = SM), 
           list(id, time), 
           env = list(id = id, time = time, SM = S * M[y])]
    n <- nrow(d)
    set(d, j = "state", value = rep_len(states, n))
    set(d, j = y, value = rep_len(rep(symbols[[y]], each = S), n))
    if (!iv(model$X_B[[y]])) {
      B <- get_B_all(
        model$gammas$gamma_B[[i]], model$X_B[[y]][1], tv(model$X_B[[y]])
      )
      d[, probability := rep_len(unlist(B), nrow(d))]
    } else {
      B <- get_B_all(
        model$gammas$gamma_B[[i]], model$X_B[[y]][ids], tv(model$X_B[[y]])
      )
      if (inherits(model, "fanhmm") && !identical(model$prior_obs, 0L)) {
        for (j in seq_along(B)) {
          B[[j]][, , 1] <- B1[[i]][, , id_idx[j]]
        }
      }
      set(d, j = "probability", value = unlist(B))
    }
    setorderv(d, c(id, time, "state"))
    out[[y]] <- d
  }
  out
}
#' @rdname emission_probs
#' @export
get_emission_probs.mnhmm <- function(model, ids = NULL, ...) {
  x <- lapply(
    split_mnhmm(model), get_emission_probs, ids = ids
  )
  stats::setNames(
    lapply(model$responses, \(y) {
      do.call(rbind, lapply(seq_along(x), \(i) {
        cbind(cluster = names(x)[i], x[[i]][[y]])
      }))
    }),
    model$responses
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
#' Prior cluster probability refers to the probability of sequence \eqn{i}
#' belonging to cluster \eqn{k} without conditioning on the observed sequence,
#' i.e., \eqn{P(i \in M_k | X_i)}, where \eqn{M_k} is cluster \eqn{k} and 
#' \eqn{X_i} are the covariates for individual \eqn{i}. The posterior cluster
#' probability refers to the probability of sequence \eqn{i} belonging to
#' cluster \eqn{k} given the observed sequence, i.e., 
#' \eqn{P(i \in M_k | y_i, X_i)}, where \eqn{y_i} is the observed sequence for 
#' individual \eqn{i}.
#' 
#' @inheritParams get_initial_probs.nhmm
#' @rdname cluster_probs
#' @param type Character string indicating the type of cluster probabilities. 
#' Either `"prior"` (the default), or `"posterior"`. See details.
#' @seealso [posterior_cluster_probabilities()].
#' @export
get_cluster_probs.mnhmm <- function(model, type = "prior", ids = NULL, ...) {
  stopifnot_(
    checkmate::test_choice(type, choices = c("prior", "posterior")), 
    "Argument {.arg type} must be either 'prior' or 'posterior'."
  )
  if (type == "posterior") {
    return(posterior_cluster_probabilities(model)[ids, ])
  }
  if (!is.null(ids)) {
    missing_ids <- setdiff(ids, unique(model$data[[model$id_variable]]))
    stopifnot_(
      length(missing_ids) == 0L,
      "The following IDs are not present in the model data: 
      {paste(missing_ids, collapse = ', ')}."
    )
    ids <- factor(ids, levels = levels(model$data[[model$id_variable]]))
  } else {
    ids <- unique(model$data[[model$id]])
  }
  if (io(model$X_omega)) {
    X <- model$X_omega[, 1L, drop = FALSE]
  } else {
    X <- model$X_omega[, ids]
  }
  D <- model$n_clusters
  N <- length(ids)
  d <- data.table(
    id = rep(ids, each = D),
    cluster = model$cluster_names,
    probability = rep_len(c(get_omega_all(model$gammas$gamma_omega, X)), D * N),
    key = "id"
  )
  setnames(d, "id", model$id_variable)
  d[]
}
#' @rdname cluster_probs
#' @export
get_cluster_probs.mhmm <- function(model, type = "prior", ...) {
  stopifnot_(
    checkmate::test_choice(type, choices = c("prior", "posterior")), 
    "Argument {.arg type} must be either 'prior' or 'posterior'."
  )
  if (type == "posterior") {
    return(posterior_cluster_probabilities(model))
  }
  pr <- exp(model$X %*% model$coefficients)
  prior_cluster_probabilities <- pr / rowSums(pr)
  if (model$n_channels == 1L) {
    ids <- rownames(model$observations)
  } else {
    ids <- rownames(model$observations[[1]])
  }
  data.table(
    cluster = model$cluster_names,
    id = rep(as_factor(ids), each = model$n_clusters),
    probability = c(t(prior_cluster_probabilities)),
    key = "id"
  )
}
