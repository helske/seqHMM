#' Log-likelihood of the Mixture Hidden Markov Model
#'
#' Function \code{logLik.mixHMModel} computes the log-likelihood value of a mixture hidden Markov model.
#'
#'
#' @export
#' @param object Hidden Markov model of class \code{mixHMModel}.
#' @param ... Ignored.
#' @return Log-likelihood of hidden Markov model.
#' @seealso \code{\link{buildMixHMM}} and \code{\link{fitMixHMM}} for building and 
#'   fitting mixture Hidden Markov models.
logLik.mixHMModel<-function(object,...){
  
  object <- combineModels(object)
  
  if(object$numberOfChannels==1){
    obsArray<-data.matrix(object$observations)-1
    obsArray[obsArray>object$numberOfSymbols]<-object$numberOfSymbols
    storage.mode(obsArray)<-"integer"
    
    logLikMixHMM(object$transitionMatrix, cbind(object$emissionMatrix,1), object$initialProbs, obsArray,
                 object$beta, object$X, object$numberOfStatesInClusters)
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
    
    logLikMixMCHMM(object$transitionMatrix, emissionArray, object$initialProbs, obsArray,
                   object$beta, object$X, object$numberOfStatesInClusters) 
    
  }
}