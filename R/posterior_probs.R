#' Posterior Probabilities for Hidden Markov Model
#'
#' Function \code{posterior_probs} computes the posterior probabilities of hidden states of
#' a hidden Markov model given the observations in logarithm scale.
#'
#' @export 
#' @param model A (mixture) hidden Markov model of class \code{hmm} or \code{mhmm}.
#' @return Posterior probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
posterior_probs<-function(model){
  
  fw <- forwardProbs(model)
  bw <- backwardProbs(model)
  ll <- logLik(model, partials = TRUE)
  
  out <- fw + bw - array(rep(ll, each = 
      sum(model$number_of_states)*model$length_of_sequences), c(sum(model$number_of_states), model$length_of_sequences, model$number_of_sequences))
  
  dimnames(out)<-list("state" = rownames(fw), "time" = 1:model$length_of_sequences, "sequence" = 1:model$number_of_sequences)
  out
}
