#' Log-likelihood of the Mixture Hidden Markov Model
#'
#' Function \code{logLik.mhmm} computes the log-likelihood value of a mixture hidden Markov model.
#'
#'
#' @export
#' @param object Hidden Markov model of class \code{mhmm}.
#' @param partials Return a vector containing the individual contributions of each sequence to the total log-likelihood. 
#' Default is FALSE, which returns the sum of all log-likelihood components.
#' @param ... Ignored.
#' @return Log-likelihood of hidden Markov model.
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_mhmm}} for building and 
#'   fitting mixture Hidden Markov models.
logLik.mhmm<-function(object, partials = FALSE, ...){
  
  object <- combine_models(object)
  
  if(object$number_of_channels == 1){
    object$observations <- list(object$observations)
    object$emission_matrix <- list(object$emission_matrix)
  }
  
  
  obsArray<-array(0,c(object$number_of_sequences,object$length_of_sequences,object$number_of_channels))
  for(i in 1:object$number_of_channels){
    obsArray[,,i]<-data.matrix(object$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>object$number_of_symbols[i]]<-object$number_of_symbols[i]
  }       
  storage.mode(obsArray)<-"integer"
  
  emissionArray<-array(1,c(object$number_of_states,max(object$number_of_symbols)+1,object$number_of_channels))
  for(i in 1:object$number_of_channels)
    emissionArray[,1:object$number_of_symbols[i],i]<-object$emission_matrix[[i]]
  
  ll <- logLikMixHMM(object$transition_matrix, emissionArray, object$initial_probs, obsArray,
    object$beta, object$X, object$number_of_states_in_clusters) 
  
  
  if(partials){
    ll
  } else sum(ll) 
}