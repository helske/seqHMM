#' Summary method for non-homogeneous hidden Markov models
#'
#' Function `summary.nhmm` gives a summary of a non-homogeneous hidden Markov model.
#'
#' @export
#' @method summary nhmm
#' @param object Non-homogeneous hidden Markov model of class `nhmm`.
summary.nhmm <- function(object, ...) {
    summary_nhmm <- list(
      transition_probs = object$transition_probs,
      emission_probs = object$emission_probs,
      initial_probs = object$initial_probs,
      logLik = ll, BIC = BIC(ll), most_probable_cluster = most_probable_cluster,
      coefficients = object$coefficients, vcov = vcov(object, conditional_se, log_space = log_space, ...),
      prior_cluster_probabilities = prior_cluster_probabilities,
      posterior_cluster_probabilities = posterior_cluster_probabilities,
      classification_table = clProbs
    )
  
  class(summary_nhmm) <- "summary.nhmm"
  summary_nhmm
}
