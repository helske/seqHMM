#' Log-likelihood of the Mixture Hidden Markov Model
#'
#' Function \code{logLik.mhmm} computes the log-likelihood value of a mixture hidden Markov model.
#'
#'
#' @export
#' @param object A mixture hidden Markov model of class \code{mhmm}.
#' @param partials Return a vector containing the individual contributions of each sequence to the total log-likelihood. 
#'   The default is \code{FALSE}, which returns the sum of all log-likelihood components.
#' @param threads Number of threads to use in parallel computing. The default is 1.
#' @param log_space Make computations using log-space instead of scaling for greater 
#' numerical stability at the cost of decreased computational performance. 
#'   The default is \code{TRUE}.
#' @param ... Ignored.
#' @return Log-likelihood of the mixture hidden Markov model. This is an object of class 
#' \code{logLik} with attributes \code{nobs} and \code{df} inherited from the model object.
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_model}} for building and 
#'   fitting mixture Hidden Markov models.
logLik.mhmm<-  function(object, partials = FALSE, threads = 1, log_space = FALSE, ...){
  
  if (threads < 1) stop("Argument threads must be a positive integer.")
  
  df <- attr(object, "df")
  nobs <- attr(object, "nobs")
  
  object <- combine_models(object)
  
  if(object$n_channels == 1){
    object$observations <- list(object$observations)
    object$emission_probs <- list(object$emission_probs)
  }
  
  
  obsArray<-array(0,c(object$n_sequences,object$length_of_sequences,object$n_channels))
  for(i in 1:object$n_channels){
    obsArray[,,i]<-data.matrix(object$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>object$n_symbols[i]]<-object$n_symbols[i]
  }       
  obsArray <- aperm(obsArray)
  
  emissionArray<-array(1,c(object$n_states,max(object$n_symbols)+1,object$n_channels))
  for(i in 1:object$n_channels)
    emissionArray[,1:object$n_symbols[i],i]<-object$emission_probs[[i]]
  
  if (!log_space) {
    ll <- logLikMixHMM(object$transition_probs, emissionArray, object$initial_probs, obsArray,
      object$coefficients, object$X, object$n_states_in_clusters, threads) 
  } else {
    ll <- log_logLikMixHMM(object$transition_probs, emissionArray, object$initial_probs, obsArray,
      object$coefficients, object$X, object$n_states_in_clusters, threads) 
  }

  
  
  structure(if (partials) ll else sum(ll), class = "logLik", df = df, nobs = nobs)
}