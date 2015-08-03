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

  if(model$numberOfChannels==1){
    obsArray<-data.matrix(model$observations)-1
    obsArray[obsArray>model$numberOfSymbols]<-model$numberOfSymbols
    storage.mode(obsArray)<-"integer"
    fw <- forward(model$transitionMatrix, cbind(model$emissionMatrix,1), 
      model$initialProbs, obsArray)    
    bw <- backward(model$transitionMatrix, cbind(model$emissionMatrix,1), obsArray)
    fw + bw
  } else{
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    for(i in 1:model$numberOfChannels){
      obsArray[,,i]<-data.matrix(model$observations[[i]])-1
      obsArray[,,i][obsArray[,,i]>model$numberOfSymbols[i]]<-model$numberOfSymbols[i]
    }    
    storage.mode(obsArray)<-"integer"
    emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    out<-forwardMC(model$transitionMatrix, emissionArray, 
      model$initialProbs, obsArray)
  }
  dimnames(out)<-list("state" = model$stateNames,"time" = 1:model$lengthOfSequences)
  out
}
