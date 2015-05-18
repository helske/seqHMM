# Divide a mixture HMM to a list of separate HMMs (covariates ignored)

divideModels <- function(model){
  
  divmodels <- replicate(model$numberOfModels, list())
  
  for(i in 1:model$numberOfModels){
    divmodels[[i]] <- buildHMM(observations=model$observations,
                               transitionMatrix=model$transitionMatrix[[i]],
                               emissionMatrix=model$emissionMatrix[[i]],
                               initialProbs=model$initialProbs[[i]])
  }
  divmodels
}
