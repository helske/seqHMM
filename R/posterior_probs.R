#' Posterior Probabilities for Hidden Markov Model
#'
#' Function \code{posterior_probs} computes the posterior probabilities of hidden states of
#' a hidden Markov model.
#'
#' @export 
#' @param model A (mixture) hidden Markov model of class \code{hmm} or \code{mhmm}.
#' @return Posterior probabilities. In case of multiple observations,
#' these are computed independently for each sequence.
posterior_probs <- function(model){
  fb <- forward_backward(model)
  fb$forward_probs * fb$backward_probs
}
