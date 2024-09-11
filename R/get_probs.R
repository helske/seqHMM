#' Get the Estimated Initial, Transition, and Emission Probabilities for NHMM 
#' or MNHMM
#' 
#' @param model An object of class `nhmm` or `mnhmm`.
#' @param newdata An optional data frame containing the new data to be used in 
#' computing the probabilities.
#' @param nsim Non-negative integer defining the number of samples from the 
#' normal approximation of the model parameters used in 
#' computing the approximate quantiles of the estimates. If `0`, only point 
#' estimates are returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`.
#' @param ... Ignored.
#' @rdname get_probs
#' @export
get_probs <- function(model, ...) {
  UseMethod("get_probs", model)
}
#' @rdname get_probs
#' @export
get_probs.nhmm <- function(model, newdata = NULL, nsim = 0, 
                           probs = c(0.025, 0.5, 0.975), ...) {
  out <- predict(model, newdata, nsim, probs, return_samples = FALSE)
  
  if (nsim > 0) {
    for(i in seq_along(probs)) {
      q <- paste0("q", 100 * probs[i])
      out$initial_probs[q] <- out$quantiles$quantiles_pi[, i]
      out$transition_probs[q] <- out$quantiles$quantiles_A[, i]
      out$emission_probs[q] <- out$quantiles$quantiles_B[, i]
    }
  }
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
get_probs.mnhmm <- function(model, newdata = NULL, nsim = 0, 
                            probs = c(0.025, 0.5, 0.975), ...) {
  
  out <- predict(model, newdata, nsim, probs, return_samples = FALSE)
  
  if (nsim > 0) {
    for(i in seq_along(probs)) {
      q <- paste0("q", 100 * probs[i])
      out$initial_probs[[q]] <- out$quantiles$quantiles_pi[, i]
      out$transition_probs[[q]] <- out$quantiles$quantiles_A[, i]
      out$emission_probs[[q]] <- out$quantiles$quantiles_B[, i]
      out$cluster_probs[[q]] <- out$quantiles$quantiles_omega[, i]
    }
  }
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
