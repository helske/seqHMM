#' Get State Names of Hidden Markov  Model
#'
#' @param object An object of class `hmm`, `mhmm`, `nhmm`, or `mnhmm`.
#' @return A character vector containing the state names, or a list of such
#'   vectors in case of mixture models.
#' @rdname state_names
#' @export
state_names <- function(object) {
  UseMethod("state_names")
}
#' @rdname state_names
#' @export
state_names.hmm <- function(object) {
  object$state_names
}
#' @rdname state_names
#' @export
state_names.mhmm <- function(object) {
  object$state_names
}
#' @rdname state_names
#' @export
state_names.nhmm <- function(object) {
  object$state_names
}
#' @rdname state_names
#' @export
state_names.mnhmm <- function(object) {
  object$state_names
}
#' Set State Names of Hidden Markov  Model
#'
#' @param object object An object of class `hmm`, `mhmm`, `nhmm`, or `mnhmm`.
#' @param value A character vector containing the new state names, or a list of
#'   such vectors in case of mixture models.
#' @return The original object with updated state names.
#' @rdname state_names
#' @export
`state_names<-` <- function(object, value) {
  UseMethod("state_names<-")
}
#' @rdname state_names
#' @export
`state_names<-.hmm` <- function(object, value) {
  stopifnot_(
    length(value) == object$n_states,
    "Number of state names does not match with the number of states."
  )
  value <- as_factor(value)
  object$state_names <- value
  names(object$initial_probs) <- value
  dimnames(object$transition_probs) <- list(from = value, to = value)
  if (object$n_channels > 1) {
    for (i in 1:object$n_channels) {
      dimnames(object$emission_probs[[i]])$state_names <- value
    }
  } else {
    dimnames(object$emission_probs)$state_names <- value
  }
  object
}
#' @rdname state_names
#' @export
`state_names<-.mhmm` <- function(object, value) {
  stopifnot_(
    length(value) == object$n_clusters,
    "New state names should be a {.cls list} with length of 
    {object$n_clusters}."
  )
  for (i in seq_len(object$n_clusters)) {
    stopifnot_(
      length(value[[i]]) == object$n_states[i],
      "Number of new state names for cluster {i} is not equal to the number of 
      hidden states."
    )
    value[[i]] <- as_factor(value[[i]])
    object$state_names[[i]] <- value[[i]]
    names(object$initial_probs[[i]]) <- value[[i]]
    dimnames(object$transition_probs[[i]]) <- list(from = value[[i]], to = value[[i]])
    if (object$n_channels > 1) {
      for (j in 1:object$n_channels) {
        dimnames(object$emission_probs[[i]][[j]])$state_names <- value[[i]]
      }
    } else {
      dimnames(object$emission_probs[[i]])$state_names <- value[[i]]
    }
  }
  object
}
#' @rdname state_names
#' @export
`state_names<-.nhmm` <- function(object, value) {
  stopifnot_(
    length(value) == object$n_states,
    "Number of state names does not match with the number of states."
  )
  object$state_names <- as_factor(value)
  object
}
#' @rdname state_names
#' @export
`state_names<-.mnhmm` <- function(object, value) {
  stopifnot_(
    length(value) == object$n_clusters,
    "New state names should be a {.cls list} with length of 
    {object$n_clusters}."
  )
  for (i in seq_len(object$n_clusters)) {
    stopifnot_(
      length(value[[i]]) == object$n_states,
      "Number of new state names for cluster {i} is not equal to the number of 
      hidden states."
    )
    object$state_names[[i]] <- as_factor(value[[i]])
  }
  object
}