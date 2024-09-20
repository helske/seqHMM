#' Get the Estimated Initial, Transition, and Emission Probabilities for NHMM 
#' or MNHMM
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param newdata An optional data frame containing the new data to be used in 
#' computing the probabilities.
#' @param ... Ignored.
#' @rdname get_probs
#' @export
get_probs <- function(model, ...) {
  UseMethod("get_probs", model)
}
#' @rdname get_probs
#' @export
get_probs.nhmm <- function(model, newdata = NULL, ...) {
  out <- predict(model, newdata)
  rownames(out$initial_probs) <- NULL
  rownames(out$transition_probs) <- NULL
  rownames(out$emission_probs) <- NULL
  list(
    initial_probs = out$initial_probs, 
    transition_probs = remove_voids(model, out$transition_probs),
    emission_probs = remove_voids(model, out$emission_probs)
  )
}
#' @rdname get_probs
#' @export
get_probs.mnhmm <- function(model, newdata = NULL, ...) {
  
  out <- predict(model, newdata)
  rownames(out$initial_probs) <- NULL
  rownames(out$transition_probs) <- NULL
  rownames(out$emission_probs) <- NULL
  rownames(out$cluster_probs) <- NULL
  list(
    initial_probs = out$initial_probs, 
    transition_probs = remove_voids(model, out$transition_probs),
    emission_probs = remove_voids(model, out$emission_probs),
    cluster_probs = out$cluster_probs
  )
}
