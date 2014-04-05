#' Most Probable Path of Hidden States of Hidden Markov Model given the sequence.
#'
#' Function \code{mostProbablePath} computes the most probable path of the hidden states of the 
#' hidden Markov model given the observed sequence.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @param tol Tolerance parameter defining the smallest possible nonzero probability. 
#' Used for ensuring that there won't be infinite log-probabilities. Default is 1e-15.
#' @return Most probable path of states given the observations. In case of multiple observations,
#' most probable path is computed independently for each sequence.
mostProbablePath<-function(model,tol=1e-15){
  model$transitionMatrix[model$transitionMatrix<tol]<-.Machine$double.xmin
  model$transitionMatrix <- log(model$transitionMatrix)
  if(is.list(model$emissionMatrix)){
    model$emissionMatrix<-lapply(model$emissionMatrix,function(x){
      x[x<tol]<-.Machine$double.xmin
      log(x)
    })
    
  } else {
    model$emissionMatrix[model$emissionMatrix<tol]<-.Machine$double.xmin
    model$emissionMatrix <- log(model$emissionMatrix)
  }
  model$initialProbs[model$initialProbs<tol]<-.Machine$double.xmin
  model$initialProbs <- log(model$initialProbs)
  
  
  if(model$numberOfChannels==1){
    
    if(model$numberOfSequences==1){
      obs<-as.integer(as.factor(model$observations))
      miss<-is.na(obs)
      #if(miss[1]) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("viterbi",PACKAGE="MVHMM",NAOK = TRUE,
                    model$transitionMatrix,model$emissionMatrix,model$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,model$lengthOfSequences,miss
                    ,q=integer(model$lengthOfSequences))
      mpp<-rownames(model$transitionMatrix)[out$q]
    } else {
      obs<-apply(model$observations,2,factor,labels=1:model$numberOfSymbols,
                 levels=model$symbolNames)
      storage.mode(obs)<-"integer"
      miss<-is.na(obs)
      #if(any(miss[,1])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      q<-matrix(0,model$numberOfSequences,model$lengthOfSequences)
      storage.mode(q)<-"integer"
      out<-.Fortran("mvviterbi",PACKAGE="MVHMM",NAOK = TRUE,
                    model$transitionMatrix,model$emissionMatrix,model$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,
                    model$lengthOfSequences,miss,model$numberOfSequences,
                    q=q)
      mpp<-apply(out$q,2,function(x) rownames(model$transitionMatrix)[x])
    }
  } else {
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    maxNumberOfSymbols<-max(model$numberOfSymbols)
    emissionArray<-array(NA,c(model$numberOfStates,maxNumberOfSymbols,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    if(model$numberOfSequences==1){
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-as.integer(as.factor(model$observations[[i]]))
      miss<-is.na(obsArray)
      #if(any(miss[1,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mcviterbi",PACKAGE="MVHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,model$lengthOfSequences,miss
                    ,q=integer(model$lengthOfSequences),model$numberOfChannels)
      mpp<-rownames(model$transitionMatrix)[out$q]
    } else {
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-apply(model$observations[[i]],2,factor,labels=1:model$numberOfSymbols[[i]],
                             levels=model$symbolNames[[i]])
      storage.mode(obsArray)<-"integer"
      miss<-is.na(obsArray)
      #if(any(miss[,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      q<-matrix(0,model$numberOfSequences,model$lengthOfSequences)
      storage.mode(q)<-"integer"
      out<-.Fortran("mvmcviterbi",PACKAGE="MVHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,
                    model$lengthOfSequences,miss,model$numberOfSequences,
                    q=q,model$numberOfChannels)
      mpp<-apply(out$q,2,function(x) rownames(model$transitionMatrix)[x])
    }
    
  }
  mpp
}
