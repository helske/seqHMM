#' Get state names from hmm or mhmm object
#'
#' @param object An object of class `hmm` or `mhmm`.
#' @return A character vector containing the state names, or a list of such
#'   vectors in `mhmm` case.
#' @export
state_names <- function(object) {
  UseMethod("state_names")
}

#' @export
state_names.hmm <- function(object) {
  object$state_names
}

#' @export
state_names.mhmm <- function(object) {
  object$state_names
}

#' Set state names for hmm or mhmm object
#'
#' @param object An object of class `hmm` or `mhmm`.
#' @param value A character vector containing the new state names, or a list of
#'   such vectors in `mhmm` case.
#' @return The modified object with updated state names.
#' @export
`state_names<-` <- function(object, value) {
  UseMethod("state_names<-")
}

#' @export
`state_names<-.hmm` <- function(object, value) {
  stopifnot_(
    length(value) == object$n_states,
    "Number of state names does not match with the number of states."
  )
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

#' @export
`state_names<-.mhmm` <- function(object, value) {
  stopifnot_(
    length(value) != object$n_clusters,
    "New state names should be a {.cls list} with length of 
    {object$n_clusters}."
  )
  for (i in 1:object$n_clusters) {
    stopifnot_(
      length(value[[i]]) == object$n_states[i],
      "Number of new state names for cluster {i} is not equal to the number of 
      hidden states."
    )
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
