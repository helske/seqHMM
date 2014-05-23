#' Estimate parameters of Hidden Markov Model
#'
#' Function \code{fitHMM} estimates the initial state, transition and emission probabilities of 
#' hidden Markov model using numerical maximization of log-likelihood. Initial values for estimation 
#' are taken from the corresponding components of model with preservation of original zero probabilities.
#' 
#' @export 
#' @param model Hidden Markov model of class HMModel.
#' @param ... Further arguments to nlm, such as arguments controlling the 
#' amount printing and number of iterations.
#' @return List with two components: First component is output from nlm and second is the estimated model.
fitHMM<-function(model,method="nlm",criterion=c("ML","L1"),...){
  
  
  
  maxIP<-which.max(model$initialProbs)
  maxIPvalue<-model$initialProbs[maxIP]
  paramIP<-{function(x) which(x>0 & x<max(x))}(model$initialProbs)
  npIP<-length(paramIP)
  
  transNZ<-which(model$transitionMatrix>0)  
  maxTM<-(max.col(model$transitionMatrix)-1)*model$numberOfStates+1:model$numberOfStates
  maxTMvalue<-apply(model$transitionMatrix,1,max)
  paramTM<-setdiff(transNZ,maxTM)
  npTM<-length(paramTM)
  
  if(model$numberOfChannels==1){
    
    emissNZ<-which(model$emissionMatrix>0)  
    
    maxEM<-(max.col(model$emissionMatrix)-1)*model$numberOfStates+1:model$numberOfStates
    maxEMvalue<-apply(model$emissionMatrix,1,max)
    paramEM<-setdiff(emissNZ,maxEM)
    npEM<-length(paramEM)
    
    
    
    initialvalues<-c(
      if(npTM>0) model$transitionMatrix[paramTM],
      if(npEM>0) model$emissionMatrix[paramEM],
      if(npIP>0) model$initialProbs[paramIP]
    )
    initialvalues <- log(initialvalues)
    
    
    likfn<-function(pars,model,estimate=TRUE){
      pars<-exp(pars) 
      if(npTM>0){
        A<-matrix(0,model$numberOfStates,model$numberOfStates)
        A[maxTM]<-maxTMvalue     
        A[paramTM]<-pars[1:npTM]
        A<-A/rowSums(A)
        model$transitionMatrix[]<-A
      }
      if(npEM>0){
        B<-matrix(0,model$numberOfStates,model$numberOfSymbols)
        B[maxEM]<-maxEMvalue     
        B[paramEM]<-pars[(npTM+1):(npTM+npEM)]
        B<-B/rowSums(B)
        model$emissionMatrix[]<-B
      }
      
      if(npIP>0){
        model$initialProbs[maxIP]<-maxIPvalue
        model$initialProbs[paramIP]<-pars[(npTM+npEM+1):length(pars)]
        model$initialProbs[]<-model$initialProbs/sum(model$initialProbs)
      }
      
      if(estimate){
        -logLik(model)
      } else model
    }
  } else {
    
    emissNZ<-lapply(model$emissionMatrix,function(x) which(x>0))
    maxEM<-lapply(model$emissionMatrix,function(x) (max.col(x)-1)*model$numberOfStates+1:model$numberOfStates)
    
    maxEMvalue<-lapply(1:model$numberOfChannels, function(x) 
      apply(model$emissionMatrix[[x]],1,max))
    
    paramEM<-lapply(1:model$numberOfChannels,function(x) setdiff(emissNZ[[x]],maxEM[[x]]))
    npEM<-sapply(paramEM,length)
    
    
    initialvalues<-c(
      if(npTM>0) model$transitionMatrix[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$numberOfChannels,
                                    function(x) model$emissionMatrix[[x]][paramEM[[x]]])),
      if(npIP>0) model$initialProbs[paramIP]
    )
    initialvalues<-log(initialvalues)
    
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    maxNumberOfSymbols<-max(model$numberOfSymbols)
    emissionArray<-array(NA,c(model$numberOfStates,maxNumberOfSymbols,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    if(model$numberOfSequences==1){
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-as.integer(as.factor(model$observations[[i]]))
    } else {
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-apply(model$observations[[i]],2,factor,labels=1:model$numberOfSymbols[[i]],
                             levels=model$symbolNames[[i]])
    }
    miss<-is.na(obsArray)
    storage.mode(miss)<-storage.mode(obsArray)<-"integer"
    
    
    likfn<-function(pars,model,estimate=TRUE){
      pars<-exp(pars)
      if(any(!is.finite(pars)) && estimate)
        return(.Machine$double.xmax)
      if(npTM>0){
        A<-matrix(0,model$numberOfStates,model$numberOfStates)
        A[maxTM]<-maxTMvalue
        A[paramTM]<-pars[1:npTM]
        A<-A/rowSums(A)
        model$transitionMatrix[]<-A
      }
      if(sum(npEM)>0){    
        emissionArray[,1:model$numberOfSymbols[1],1][maxEM[[1]]]<-maxEMvalue[[1]]     
        emissionArray[,1:model$numberOfSymbols[1],1][paramEM[[1]]]<-pars[(npTM+1):(npTM+cumsum(npEM)[1])]
        emissionArray[,1:model$numberOfSymbols[1],1]<-
          emissionArray[,1:model$numberOfSymbols[1],1]/rowSums(emissionArray[,1:model$numberOfSymbols[1],1])
        for(i in 2:model$numberOfChannels){
          emissionArray[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]  
          emissionArray[,1:model$numberOfSymbols[i],i][paramEM[[i]]]<-pars[(npTM+1+cumsum(npEM)[i-1]):(npTM+cumsum(npEM)[i])]
          emissionArray[,1:model$numberOfSymbols[i],i]<-
            emissionArray[,1:model$numberOfSymbols[i],i]/rowSums(emissionArray[,1:model$numberOfSymbols[i],i])
        }
      }
      
      if(npIP>0){
        model$initialProbs[maxIP]<-maxIPvalue
        model$initialProbs[paramIP]<-pars[(npTM+sum(npEM)+1):length(pars)]
        model$initialProbs[]<-model$initialProbs/sum(model$initialProbs)
      }
      
      if(estimate){
        if(criterion=="ML"){
          if(model$numberOfSequences==1){
            ll<-.Fortran("mchmmloglik",PACKAGE="seqHMM",NAOK = TRUE,
                         model$transitionMatrix,emissionArray,model$initialProbs,
                         obsArray,model$numberOfStates,maxNumberOfSymbols,
                         model$lengthOfSequences,miss,logLik=double(1),model$numberOfChannels)$logLik
          } else{
            ll<-.Fortran("mvmchmmloglik",PACKAGE="seqHMM",NAOK = TRUE,
                         model$transitionMatrix,emissionArray,model$initialProbs,
                         obsArray,model$numberOfStates,maxNumberOfSymbols,
                         model$lengthOfSequences,miss,model$numberOfSequences,logLik=double(1),
                         model$numberOfChannels)$logLik
          }
          -ll
        } else {
          
          n<-model$lengthOfSequences
          k<-model$numberOfSequences
          
          #polkujen todennäköisyydet (yksi jokaiselle sekvenssille)
          viterbi<-mostProbablePath(model)$mpp
          
          #pStates<-viterbi$logP
          
          #sekvenssien ja piilotilojen todennäköisyydet ehdolla polut
          pSeq<-pStates<-numeric(k)               
          names(model$initialProbs)<-rownames(model$transitionMatrix)
          for(j in 1:k){
            pStates[j]<-log(model$initialProbs[viterbi[j,1]])+
              sum(log(model$trans[cbind(viterbi[j,-n],viterbi[j,-1])]))
            
            for(s in 1:model$numberOfChannels)
              pSeq[j]<-pSeq[j]+sum(log(
                emissionArray[cbind(viterbi[j,],obsArray[j,,s],s)]),na.rm=TRUE)
            
          }
          
          sum(pStates)+sum(pSeq)
          
        }
      } else {
        if(sum(npEM)>0){    
          model$emissionMatrix[[1]][]<-emissionArray[,1:model$numberOfSymbols[1],1]
          for(i in 2:model$numberOfChannels){
            model$emissionMatrix[[i]][]<-emissionArray[,1:model$numberOfSymbols[i],i]
          }
        }
        model
      }
      
    }   
  }
  browser()
  #if(method=="nlm"){
  fit<-nlm(f=likfn,p=initialvalues,
           model=model,...)
  list(opt=fit,model=likfn(fit$e,model,FALSE))
  #} else{
  #  fit<-bobyqa(fn=likfn, par=initialvalues, lower=0, model=model,...)
  #  if(fit$ierr!=0)
  #    print(fit$msg)
  #  out<-list(opt=fit,model=likfn(fit$p,model,FALSE))
  #}
  #out
}