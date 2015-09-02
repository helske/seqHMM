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
  
  if(inherits(model,"mixHMModel")){
    fw <- forwardProbs(model)[,model$lengthOfSequences,]
    model <- combineModels(model)
    mix<-TRUE
  } else mix <- FALSE
  
  if(model$numberOfChannels == 1){
    model$observations <- list(model$observations)
    model$emissionMatrix <- list(model$emissionMatrix)
  }
  
  obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
  for(i in 1:model$numberOfChannels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$numberOfSymbols[i]]<-model$numberOfSymbols[i]
  }    
  storage.mode(obsArray)<-"integer"
  emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
  for(i in 1:model$numberOfChannels)
    emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
  
  fw <- forward(model$transitionMatrix, emissionArray, 
    model$initialProbs, obsArray)
  bw <- backward(model$transitionMatrix, emissionArray, obsArray)
  ll <- logLikHMM(model$transitionMatrix, emissionArray, 
    model$initialProbs, obsArray)
  out <- fw + bw - ll
  
  dimnames(out)<-list("state" = model$stateNames,"time" = 1:model$lengthOfSequences, "sequence" = 1:model$numberOfSequences)
  out
}
