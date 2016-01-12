#' Posterior Probabilities for (Mixture) Hidden Markov Models
#'
#' Function \code{posterior_probs} computes the posterior probabilities of hidden states of
#' a (mixture) hidden Markov model.
#'
#' @export 
#' @param model A (mixture) hidden Markov model of class \code{hmm} or \code{mhmm}.
#' @param log_space Compute posterior probabilities in logarithmic scale. The default is \code{FALSE}.
#' @return Posterior probabilities. In case of multiple observations,
#' these are computed independentlsy for each sequence.
#' @examples 
#' # Load a pre-defined MHMM
#' data("mhmm_biofam")
#' 
#' # Compute posterior probabilities
#' pb <- posterior_probs(mhmm_biofam)
#' 
#' # Locally most probable states for the first subject:
#' pb[, , 1]
posterior_probs <- function(model, log_space = FALSE){
  fb <- forward_backward(model, log_space = log_space)
  if (!log_space) {
  fb$forward_probs * fb$backward_probs
  } else {
    ll <- logLik(model, partials = TRUE, log_space = TRUE)
    fb$forward_probs + fb$backward_probs - 
      array(rep(ll, each = sum(model$n_states) * model$length_of_sequences), 
        c(sum(model$n_states), model$length_of_sequences, model$n_sequences))
  }
}
