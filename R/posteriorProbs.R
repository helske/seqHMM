#' Posterior Probabilities for Hidden Markov Model
#'
#' Function \code{posteriorProbs} computes the posterior probabilities of hidden states of
#' hidden Markov model given the observations in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @return Posterior probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
posteriorProbs<-function(model){
  
  fw <- forwardProbs(model)
  bw <- backwardProbs(model)
  ll <- logLik(model, partials = TRUE)
  
  out <- fw + bw - array(rep(ll, each = 
      sum(model$numberOfStates)*model$lengthOfSequences), c(sum(model$numberOfStates), model$lengthOfSequences, model$numberOfSequences))
  
  dimnames(out)<-list("state" = rownames(fw), "time" = 1:model$lengthOfSequences, "sequence" = 1:model$numberOfSequences)
  out
}
