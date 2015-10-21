#' Posterior Probabilities for Hidden Markov Model
#'
#' Function \code{posterior_probs} computes the posterior probabilities of hidden states of
#' a hidden Markov model.
#'
#' @export 
#' @param model A (mixture) hidden Markov model of class \code{hmm} or \code{mhmm}.
#' @param log_space Compute posterior probabilities in logarithmic scale. 
#' @return Posterior probabilities. In case of multiple observations,
#' these are computed independently for each sequence.
posterior_probs <- function(model, log_space = FALSE){
  fb <- forward_backward(model, log_space = log_space)
  if (!log_space) {
  fb$forward_probs * fb$backward_probs
  } else {
    ll <- logLik(model, partials = TRUE, log_space = TRUE)
    fb$forward_probs + fb$backward_probs - array(rep(ll, each = 
        sum(model$n_states)*model$length_of_sequences), c(sum(model$n_states), model$length_of_sequences, model$n_sequences))
  }
}
