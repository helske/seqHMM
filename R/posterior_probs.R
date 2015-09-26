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
  
  post_probs <- fb$forward_probs * fb$backward_probs
  
  if (!is.null(time_names <- colnames(model$observations[[1]]))) {
    time_names <- 1:model$length_of_sequences
  }
  
  if (!is.null(sequence_names <- rownames(model$observations[[1]]))) {
    sequence_names <- 1:model$n_sequences
  } 
  
  dimnames(post_probs) <- list("state" = model$state_names, 
    "time" = time_names, "sequence" = sequence_names)
  post_probs
}
