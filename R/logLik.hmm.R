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
  
  if(object$n_channels == 1){
    object$observations <- list(object$observations)
    object$emission_matrix <- list(object$emission_matrix)
  }
  
  
  obsArray<-array(0,c(object$n_sequences,object$length_of_sequences,object$n_channels))
  for(i in 1:object$n_channels){
    obsArray[,,i]<-data.matrix(object$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>object$n_symbols[i]]<-object$n_symbols[i]
  }       
  storage.mode(obsArray)<-"integer"
  
  emissionArray<-array(1,c(object$n_states,max(object$n_symbols)+1,object$n_channels))
  for(i in 1:object$n_channels)
    emissionArray[,1:object$n_symbols[i],i]<-object$emission_matrix[[i]]
  
  ll <- logLikHMM(object$transition_matrix, emissionArray, 
    object$initial_probs, obsArray)
  
  structure(if (partials) ll else sum(ll), class = "logLik", df = attr(object, "df"), nobs = attr(object, "nobs"))
  
}