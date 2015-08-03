#' Forward Probabilities for Hidden Markov Model
#'
#' Function \code{forwardProbs} computes the forward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @return Forward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
forwardProbs<-function(model){
  if(model$numberOfChannels==1){
    obsArray<-data.matrix(model$observations)-1
    obsArray[obsArray>model$numberOfSymbols]<-model$numberOfSymbols
    storage.mode(obsArray)<-"integer"
    out<-forward(model$transitionMatrix, cbind(model$emissionMatrix,1), 
                 model$initialProbs, obsArray)
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


#' Backward Probabilities for Hidden Markov Model
#'
#' Function \code{backwardProbs} computes the backward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @return Backward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
backwardProbs<-function(model){
  if(model$numberOfChannels==1){
    obsArray<-data.matrix(model$observations)-1
    obsArray[obsArray>model$numberOfSymbols]<-model$numberOfSymbols
    storage.mode(obsArray)<-"integer"
    out<-backward(model$transitionMatrix, cbind(model$emissionMatrix,1), obsArray)
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
    
    out<-backwardMC(model$transitionMatrix, emissionArray, obsArray)
  }
  dimnames(out)<-list("state" = model$stateNames,"time" = 1:model$lengthOfSequences)
  out
  
}
