#' Log-likelihood of the Hidden Markov Model
#'
#' Function \code{logLik.hmm} computes the log-likelihood value of a hidden Markov model.
#'
#'
#' @export
#' @importFrom stats logLik
#' @param object Hidden Markov model of class \code{hmm}.
#' @param partials Return a vector containing the individual contributions of each sequence to the total log-likelihood. 
#' Default is FALSE, which returns the sum of all log-likelihood components.
#' @param ... Ignored.
#' @return Log-likelihood of hidden Markov model.
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting Hidden Markov models.
logLik.hmm<-function(object, partials = FALSE, ...){
  
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
  
  ll <- logLikHMM(object$transition_matrix, emissionArray, 
    object$initial_probs, obsArray)
  
  if(partials){
    ll
  } else sum(ll) 
}