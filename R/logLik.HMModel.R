#' Log-likelihood of the Hidden Markov Model
#'
#' Function \code{logLik.HMModel} computes the log-likelihood value of a hidden Markov model.
#'
#'
#' @export
#' @param object Hidden Markov model of class \code{HMModel}.
#' @param ... Ignored.
#' @return Log-likelihood of hidden Markov model.
#' @seealso \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and 
#'   fitting Hidden Markov models.
logLik.HMModel<-function(object,...){
  
  if(object$numberOfChannels==1){
    obsArray<-data.matrix(object$observations)-1
    obsArray[obsArray>object$numberOfSymbols]<-object$numberOfSymbols
    storage.mode(obsArray)<-"integer"
    
      logLikHMM(object$transitionMatrix, cbind(object$emissionMatrix,1), 
                object$initialProbs, obsArray)
  } else {
    obsArray<-array(0,c(object$numberOfSequences,object$lengthOfSequences,object$numberOfChannels))
    for(i in 1:object$numberOfChannels){
      obsArray[,,i]<-data.matrix(object$observations[[i]])-1
      obsArray[,,i][obsArray[,,i]>object$numberOfSymbols[i]]<-object$numberOfSymbols[i]
    }       
    storage.mode(obsArray)<-"integer"
    
    emissionArray<-array(1,c(object$numberOfStates,max(object$numberOfSymbols)+1,object$numberOfChannels))
    for(i in 1:object$numberOfChannels)
      emissionArray[,1:object$numberOfSymbols[i],i]<-object$emissionMatrix[[i]]
    
      logLikMCHMM(object$transitionMatrix, emissionArray, 
                  object$initialProbs, obsArray)
    
  }
}