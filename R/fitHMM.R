#' Estimate parameters of Hidden Markov Model
#'
#' Function \code{fitHMM} estimates the initial state, transition and emission probabilities of 
#' hidden Markov model using numerical maximization of log-likelihood. Initial values for estimation 
#' are taken from the corresponding components of model with preservation of original zero probabilities.
#' 
#' If both maxit1 and maxit2 are positive, function first starts
#' 
#' @export 
#' @param model Hidden Markov model of class HMModel or MCHMModel.
#' @param maxit1 Number of iterations for EM algorithm.
#' @param maxit2 Number of iterations for direct numerical optimization algorithm.
#' @param ... Further arguments to \code{nlm}, such as arguments controlling the 
#' amount printing.
#' @return List with two components: First component is output from nlm and second is the estimated model.
fitHMM<-function(model,maxit1=50,maxit2=10000,reltol1=1e-8,...){
  
  if(class(model)=="HMModel"){
    obsArray<-data.matrix(model$observations)-1
    obsArray[is.na(obsArray)]<-model$numberOfSymbols
    storage.mode(obsArray)<-"integer"
  } else{
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    for(i in 1:model$numberOfChannels){
      obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    }
      obsArray[is.na(obsArray)]<-max(model$numberOfSymbols)
    
    storage.mode(obsArray)<-"integer"
    
  }
  if(maxit1>0){
    if(class(model)=="HMModel"){
      
      resEM<-EM(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, 
                obsArray, model$numberOfSymbols, maxit1, reltol1)
      if(resEM$change<0)
        warning("EM algorithm stopped due to the decreasing log-likelihood. ")
      
      model$initialProbs[]<-resEM$initialProbs
      model$transitionMatrix[]<-resEM$transitionMatrix
      model$emissionMatrix[]<-resEM$emissionMatrix[,1:model$numberOfSymbols]
      
    } else {
      
      
      
      emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
      for(i in 1:model$numberOfChannels)
        emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
      
      resEM<-EMMC(model$transitionMatrix, emissionArray, model$initialProbs, obsArray, model$numberOfSymbols, maxit1, reltol1)
      if(resEM$change<0)
        warning("EM algorithm stopped due to the decreasing log-likelihood. ")
      
      model$initialProbs[]<-resEM$initialProbs
      model$transitionMatrix[]<-resEM$transitionMatrix
      model$emissionMatrix<-lapply(seq(dim(resEM$emissionArray)[3]), 
                                   function(i) resEM$emissionArray[ , 1:model$numberOfSymbols[i], i])
      names(model$emissionMatrix)<-model$channelNames
    }
    
  } else resEM <-NULL
  if(maxit2>0){
    maxIP<-which.max(model$initialProbs)
    maxIPvalue<-model$initialProbs[maxIP]
    paramIP<-{function(x) which(x>0 & x<max(x))}(model$initialProbs)
    npIP<-length(paramIP)
    
    transNZ<-which(model$transitionMatrix>0)  
    maxTM<-(max.col(model$transitionMatrix)-1)*model$numberOfStates+1:model$numberOfStates
    maxTMvalue<-apply(model$transitionMatrix,1,max)
    paramTM<-setdiff(transNZ,maxTM)
    npTM<-length(paramTM)
    
    if(class(model)=="HMModel"){
      
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
          
        }
        if(npEM>0){
          B<-matrix(0,model$numberOfStates,model$numberOfSymbols)
          B[maxEM]<-maxEMvalue     
          B[paramEM]<-pars[(npTM+1):(npTM+npEM)]
          B<-cbind(B/rowSums(B),1)          
        }
        
        if(npIP>0){
          model$initialProbs[maxIP]<-maxIPvalue
          model$initialProbs[paramIP]<-pars[(npTM+npEM+1):length(pars)]
          model$initialProbs[]<-model$initialProbs/sum(model$initialProbs)
        }
        
        if(estimate){
          - logLikHMM(A, B, model$initialProbs, obsArray)
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

      
      emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
      for(i in 1:model$numberOfChannels)
        emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
      
      
      
      likfn<-function(pars,model,estimate=TRUE){
        pars<-exp(pars)
        if(any(!is.finite(pars)) && estimate)
          return(.Machine$double.xmax)
        if(npTM>0){
          A<-matrix(0,model$numberOfStates,model$numberOfStates)
          A[maxTM]<-maxTMvalue
          A[paramTM]<-pars[1:npTM]
          A<-A/rowSums(A)
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
          - logLikMCHMM(A, emissionArray, 
                        model$initialProbs, obsArray)
          
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
    
    
    fit<-nlm(f=likfn,p=initialvalues, model=model,iterlim=maxit2,...)
    model=likfn(fit$e,model,FALSE)
  } else fit<-NULL
  list(model=model,EM=resEM,ML=fit)
}