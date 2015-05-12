combineModels <- function(model){
  
  numberOfStatesInModels <- model$numberOfStates
  numberOfStates <- sum(model$numberOfStates)
  transitionMatrix <- as.matrix(.bdiag(model$transitionMatrix))
  stateNames <- unlist(originalStateNames<-model$stateNames)
  if (length(unique(stateNames))!= length(stateNames)){
    stateNames <- paste(stateNames,rep(1:model$numberOfModels,model$numberOfStates),sep="_")
  }
  dimnames(transitionMatrix) <- replicate(2, stateNames, simplify=FALSE)
  
  if(model$numberOfChannels>1){
    emissionMatrix <- lapply(1:model$numberOfChannels, function(i){
      x <- do.call("rbind",sapply(model$emissionMatrix,"[",i))
      rownames(x) <- stateNames
      x
    })
    names(emissionMatrix) <- model$channelNames
  } else {
    emissionMatrix <- do.call("rbind",model$emissionMatrix)
    rownames(emissionMatrix) <- stateNames
  }
  
  model <- list(observations = model$observations,
                transitionMatrix = transitionMatrix,
                emissionMatrix=emissionMatrix,
                initialProbs = unlist(model$initialProbs),
                beta=model$beta, X=model$X,
                modelNames=model$modelNames,
                stateNames=stateNames,
                symbolNames=model$symbolNames,
                channelNames=model$channelNames,
                lengthOfSequences=model$lengthOfSequences,
                numberOfSequences=model$numberOfSequences,
                numberOfSymbols=model$numberOfSymbols,
                numberOfStates=numberOfStates,
                numberOfChannels=model$numberOfChannels,
                numberOfCovariates=model$numberOfCovariates,
                numberOfModels=model$numberOfModels,
                numberOfStatesInModels=numberOfStatesInModels,
                originalStateNames = originalStateNames)
  class(model)<-"combined_mixHMModel"
  model
}