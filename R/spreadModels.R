# Transform a mixHMModel object to separate HMModel objects

spreadModels <- function(model){
  
  stnames <- substr(model$stateNames, 1, 1)
  rownames(model$transitionMatrix) <- colnames(model$transitionMatrix) <- stnames
  
  for(j in 1:model$numberOfChannels){
    rownames(model$emissionMatrix[[j]]) <- stnames
  }
  
  names(model$emissionMatrix) <- model$channelNames
  
  transM <- vector("list", model$numberOfModels)
  emissM <- vector("list", model$numberOfModels)
  init <- vector("list", model$numberOfModels)
  k <- 0
  for(m in 1:model$numberOfModels){
    transM[[m]] <- model$transitionMatrix[(k+1):(k+model$numberOfStatesInModels[m]),
                                          (k+1):(k+model$numberOfStatesInModels[m])]
    for(j in 1:model$numberOfChannels){
      emissM[[m]][[j]] <- model$emissionMatrix[[j]][(k+1):(k+model$numberOfStatesInModels[m]),]
    }
    names(emissM[[m]]) <- model$channelNames
    init[[m]] <- model$initialProbs[(k+1):(k+model$numberOfStatesInModels[m])]
    k <- sum(model$numberOfStatesInModels[1:m])
  }
  
  names(transM) <- names(emissM) <- names(init) <- model$modelNames
  
  model$transitionMatrix <- transM
  model$emissionMatrix <- emissM
  model$initialProbs <- init
  class(model) <-"mixHMModel"
  model
}