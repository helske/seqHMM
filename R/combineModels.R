combineModels <- function(model){
  
  numberOfStates <- sum(model$numberOfStates)
  transitionMatrix <- as.matrix(.bdiag(model$transitionMatrix))
  stateNames <- unlist(model$stateNames)
  if (length(unique(stateNames))!= length(stateNames)){
    stateNames <- paste(stateNames,rep(1:model$numberOfModels,model$numberOfStates),sep=".")
  }
  dimnames(transitionMatrix) <- replicate(2, stateNames, simplify=FALSE)
  
  emissionMatrix <- lapply(1:model$numberOfChannels, function(i){
    x <- do.call("rbind",sapply(model$emissionMatrix,"[",i))
    rownames(x) <- stateNames
    x
  })
  names(emissionMatrix) <- model$channelNames
  
  model <- list(observations = model$observations,
                transitionMatrix = transitionMatrix,
                emissionMatrix=emissionMatrix,
                initialProbs <- unlist(model$initialProbs),
                beta=beta, X=X,
                stateNames=stateNames,
                symbolNames=symbolNames,channelNames=model$channelNames,lengthOfSequences=model$lengthOfSequences,
                numberOfSequences=model$numberOfSequences,
                numberOfSymbols=model$numberOfSymbols,numberOfStates=numberOfStates,
                numberOfChannels=model$numberOfChannels,
                numberOfCovariates=numberOfCovariates)
  class(model)<-"mixHMModel"
  model
}