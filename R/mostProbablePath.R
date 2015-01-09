#' Most Probable Path of Hidden States of Hidden Markov Model given the sequence.
#'
#' Function \code{mostProbablePath} computes the most probable path of the hidden states of the 
#' hidden Markov model given the observed sequence.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel} or \code{MCHMModel}.
#' @return List which contains the most probable path of states (mpp) given the observations and its 
#' log-probability (logP). In case of multiple observations, most probable path is computed independently for each sequence.
mostProbablePath<-function(model){
  
  if(model$numberOfCovariates==0){
  model$initialProbs <- log(model$initialProbs)
  model$initialProbs[!is.finite(model$initialProbs)]<- -0.1*.Machine$double.xmax
  }
  model$transitionMatrix <- log(model$transitionMatrix)
  model$transitionMatrix[!is.finite(model$transitionMatrix)]<- -0.1*.Machine$double.xmax
  if(model$numberOfChannels==1){
    
    model$emissionMatrix <- log(model$emissionMatrix)
    model$emissionMatrix[!is.finite(model$emissionMatrix)]<- -0.1*.Machine$double.xmax
    obsArray<-data.matrix(model$observations)-1
    obsArray[obsArray>model$numberOfSymbols]<-model$numberOfSymbols
    storage.mode(obsArray)<-"integer"
    if(model$numberOfCovariates==0){
    out<-viterbi(model$transitionMatrix, cbind(model$emissionMatrix,0), 
                 model$initialProbs, obsArray)
    } else out<-viterbix(model$transitionMatrix, cbind(model$emissionMatrix,0), 
                        obsArray,model$beta, model$X)
    if(model$numberOfSequences==1){
      mpp<-t(rownames(model$transitionMatrix)[out$q+1])
    }else{
      mpp<-apply(out$q+1,2,function(x) rownames(model$transitionMatrix)[x])
    }
    mpp<-seqdef(mpp,alphabet=model$stateNames,
          id=rownames(model$obs),
          start=attr(model$obs,"start"),
          xtstep=attr(model$obs,"xtstep"))
  } else {
    model$emissionMatrix<-lapply(model$emissionMatrix,function(x){
      x<-log(x)
      x[!is.finite(x)]<- -0.1*.Machine$double.xmax
      x
    })
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    for(i in 1:model$numberOfChannels){
      obsArray[,,i]<-data.matrix(model$observations[[i]])-1
      obsArray[,,i][obsArray[,,i]>model$numberOfSymbols[i]]<-model$numberOfSymbols[i]
    }       
    storage.mode(obsArray)<-"integer"
    
    emissionArray<-array(0,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    if(model$numberOfCovariates==0){
    out<-viterbiMC(model$transitionMatrix, emissionArray, 
                   model$initialProbs, obsArray)
    } else out<-viterbiMCx(model$transitionMatrix, emissionArray, 
                     obsArray, model$beta, model$X)
    
    if(model$numberOfSequences==1){
      mpp<-t(rownames(model$transitionMatrix)[out$q+1])
    }else{
      mpp<-apply(out$q+1,2,function(x) rownames(model$transitionMatrix)[x])
    }
    mpp<-seqdef(mpp,alphabet=model$stateNames,
                id=rownames(model$obs[[1]]),
                start=attr(model$obs[[1]],"start"),
                xtstep=attr(model$obs[[1]],"xtstep"))
  }
  
  list(mpp=mpp,logP=out$logp)
  
}
