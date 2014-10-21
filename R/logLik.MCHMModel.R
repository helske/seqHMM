#' Log-likelihood of the Multichannel Hidden Markov Model
#'
#' Function \code{logLik.MCHMModel} computes the log-likelihood value of a multichannel hidden Markov model.
#'
#'
#' @export
#' @aliases logLik logLik.MCHMModel
#' @param object Hidden Markov model of class \code{MCHMModel}.
#' @return Log-likelihood of multichannel hidden Markov model.
logLik.MCHMModel<-function(object,...){
  
  obsArray<-array(0,c(object$numberOfSequences,object$lengthOfSequences,object$numberOfChannels))
  for(i in 1:object$numberOfChannels){
    obsArray[,,i]<-data.matrix(object$observations[[i]])-1
  }
  obsArray[is.na(obsArray)]<-max(object$numberOfSymbols)
  
  storage.mode(obsArray)<-"integer"
  
  emissionArray<-array(1,c(object$numberOfStates,max(object$numberOfSymbols)+1,object$numberOfChannels))
  for(i in 1:object$numberOfChannels)
    emissionArray[,1:object$numberOfSymbols[i],i]<-object$emissionMatrix[[i]]
  
  logLikMCHMM(object$transitionMatrix, emissionArray, 
              object$initialProbs, obsArray)
  
}