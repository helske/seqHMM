#' @export
trimHMM<-function(model,maxit=0,return.loglik=FALSE,zerotol=1e-8,...){
  
  if(!(any(model$initialProbs<zerotol & model$initialProbs>0) || 
         any(model$transitionMatrix<zerotol & model$transitionMatrix>0)
       || any(sapply(model$emissionMatrix,function(x) any(x<zerotol & x>0))))){
    print("Nothing to trim.")
    if(return.loglik){
    return(list(model=model,loglik=logLik(model)))
    } else return(model)
  }
  ll_original<- logLik(model)
  
  model$initialProbs[model$initialProbs<zerotol]<-0
  model$initialProbs<-model$initialProbs/sum(model$initialProbs)
  model$transitionMatrix[model$transitionMatrix<zerotol]<-0
  model$transitionMatrix<-model$transitionMatrix/rowSums(model$transitionMatrix)
  if(model$numberOfChannels==1){
    model$emissionMatrix[model$emissionMatrix<zerotol]<-0
    model$emissionMatrix<-model$emissionMatrix/rowSums(model$emissionMatrix)
    
  }else {
    for(i in 1:model$numberOfChannels){
      model$emissionMatrix[[i]][model$emissionMatrix[[i]]<zerotol]<-0
      model$emissionMatrix[[i]]<-model$emissionMatrix[[i]]/
        rowSums(model$emissionMatrix[[i]])
      
    }
  }
  if(!is.finite(ll0<-logLik(model)))
    stop("Initial trimming resulted a non-finite log-likelihood. Try changing zerotol parameter.")
  
  if(maxit>0){
    for(ii in 1:maxit){
      fit<-fitHMM(model,...)
      ll<- -fit$opt$m
      
      if(ll>ll0){
        model<-fit$model
        ll0<-ll
      } else break
      
      if(!(any(model$initialProbs<zerotol & model$initialProbs>0) || 
             any(model$transitionMatrix<zerotol & model$transitionMatrix>0)
         || any(sapply(model$emissionMatrix,function(x) any(x<zerotol & x>0)))))
        break
      
      model$initialProbs[model$initialProbs<zerotol]<-0
      model$initialProbs<-model$initialProbs/sum(model$initialProbs)
      model$transitionMatrix[model$transitionMatrix<zerotol]<-0
      model$transitionMatrix<-model$transitionMatrix/rowSums(model$transitionMatrix)
      if(model$numberOfChannels==1){
        model$emissionMatrix[model$emissionMatrix<zerotol]<-0
        model$emissionMatrix<-model$emissionMatrix/rowSums(model$emissionMatrix)
        
      }else {
        for(i in 1:model$numberOfChannels){
          model$emissionMatrix[[i]][model$emissionMatrix[[i]]<zerotol]<-0
          model$emissionMatrix[[i]]<-model$emissionMatrix[[i]]/
            rowSums(model$emissionMatrix[[i]])
          
        }
      }
    }
    print(ii)
  }
  if(maxit>0)
  print(paste(ii,"iterations used."))
  
  if(ll0<ll_original){
    print(paste("Log-likelihood of trimmed model is smaller than the original log-likelihood, ll_trim-ll_orig =", ll0-ll_original))
  } else print(paste("Trimming improved log-likelihood, ll_trim-ll_orig =", ll0-ll_original))
  
 
  if(return.loglik){
    list(model=model,loglik=ll0)
  } else model 
}