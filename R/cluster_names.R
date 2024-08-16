#' Get cluster names from mhmm object
#'
#' @param object An object of class `mhmm`.
#' @return A character vector containing the cluster names.
#' @export
cluster_names <- function(object) {
  UseMethod("cluster_names")
}

#' @export
cluster_names.mhmm <- function(object) {
  object$cluster_names
}

#' Set Cluster Names for Mixture Models
#'
#' @param object An object of class `mhmm` or `mnhmm`.
#' @param value A character vector containing the new cluster names.
#' @return The modified object with updated cluster names.
#' @export
`cluster_names<-` <- function(object, value) {
  UseMethod("cluster_names<-")
}

#' @export
`cluster_names<-.mnhmm` <- function(object, value) {
  if (length(value) != object$n_clusters) {
    stop(
      paste0(
        "New cluster names should be a vector of length",
        object$n_clusters, "."
      )
    )
  } else {
    object$cluster_names <- value
    names(object$state_names) <- value
  }
  object
}
#' @export
`cluster_names<-.mhmm` <- function(object, value) {
  if (length(value) != object$n_clusters) {
    stop(
      paste0(
        "New cluster names should be a vector of length",
        object$n_clusters, "."
      )
    )
  } else {
    object$cluster_names <- value
    names(object$state_names) <- value
    colnames(object$coefficients) <- value
    names(object$transition_probs) <- names(object$emission_probs) <-
      names(object$initial_probs) <- value
  }
  object
}