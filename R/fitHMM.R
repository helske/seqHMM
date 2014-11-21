#' Estimate parameters of Hidden Markov Model
#'
#' Function \code{fitHMM} estimates the initial state, transition and emission probabilities of 
#' hidden Markov model using numerical maximization of log-likelihood. Initial values for estimation 
#' are taken from the corresponding components of model with preservation of original zero probabilities.
#' 
#' By default, estimation start with EM algorithm and then switches to direct numerical maximization.
#' 
#' @export 
#' @import optimx
#' @param model Hidden Markov model of class HMModel or MCHMModel.
#' @param use.em Logical, use EM algorithm at the start of parameter estimation. Default is TRUE.
#' @param use.optimx Logical, use direct numerical optimization via \code{\link{optimx}} after EM algorithm. Default is TRUE. 
#' @param em.control Optional list of control parameters for for EM algorithm. Possible arguments are 
#' \itemize{
#'  \item{maxit}{Maximum number of iterations, default is 100.}
#'  \item{trace}{Level of printing. Possible values are 0 (prints nothing),
#'   1 (prints information at start and end of algorithm), and 2 (prints at every iteration).}
#'   \item{reltol}{Relative tolerance for convergence defined as \eqn{(tmp - sumlogLik)/(abs(sumlogLik)+0.1)}. Default is 1e-8.}
#' }
#' @param method Optimization method used by \code{optimx}. Default is \code{"BFGS"}. Note that \code{fitHMM} uses 
#' Softmax parameterization so unconstrained optimization methods are used.
#' @param itnmax Maximum number of iterations use by \code{optimx}. Default is 10000.
#' @param optimx.control Optional list of additional arguments for \code{\link{optimx}} argument \code{control}. 
#' Note that default values for \code{starttests} and \code{kkt} are set to \code{FALSE}, which differs from the default
#' behaviour of \code{optimx}. If EM algorithm is used, \code{fnscale} is also set to current optimum (unless modified by user).
#' @param ... Additional arguments to optimx.
#' @return List with components
#' \itemize{
#'  \item{model}{Estimated model.}
#'  \item{logLik}{Log-likelihood of the estimated model.}
#'   \item{em.results}{Results from EM algorithm.}
#'   \item{optimx.results}{Results from direct numerical optimization via \code{\link{optimx}}.}
#' }
#'
fitHMM<-function(model,use.em=TRUE,use.optimx=TRUE,em.control=list(),method="BFGS",itnmax=10000,optimx.control=list(),...){
  
  
  if(model$numberOfChannels==1){
    obsArray<-data.matrix(model$observations)-1
    obsArray[obsArray>model$numberOfSymbols]<-model$numberOfSymbols
    
  } else{
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    for(i in 1:model$numberOfChannels){
      obsArray[,,i]<-data.matrix(model$observations[[i]])-1
      obsArray[,,i][obsArray[,,i]>model$numberOfSymbols[i]]<-model$numberOfSymbols[i]
    }           
  }
  #storage.mode(obsArray)<-"integer"
  if(use.em){
    em.con <- list(trace = 0, maxit=100,reltol=1e-8)
    nmsC <- names(em.con)  
    em.con[(namc <- names(em.control))] <- em.control
    if (length(noNms <- namc[!namc %in% nmsC])) 
      warning("Unknown names in em.control: ", paste(noNms, collapse = ", "))
    
    if(model$numberOfChannels==1){
      
      resEM<-EM(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, 
                obsArray, model$numberOfSymbols, em.con$maxit, em.con$reltol,em.con$trace)
      if(resEM$change<0)
        warning("EM algorithm stopped due to the decreasing log-likelihood. ")
      
      model$initialProbs[]<-resEM$initialProbs
      model$transitionMatrix[]<-resEM$transitionMatrix
      model$emissionMatrix[]<-resEM$emissionMatrix[,1:model$numberOfSymbols]
      
    } else {    
      
      emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
      for(i in 1:model$numberOfChannels)
        emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
      
      resEM<-EMMC(model$transitionMatrix, emissionArray, model$initialProbs, obsArray, 
                  model$numberOfSymbols, em.con$maxit, em.con$reltol,em.con$trace)
      if(resEM$change<0)
        warning("EM algorithm stopped due to the decreasing log-likelihood. ")
      
      model$initialProbs[]<-resEM$initialProbs
      model$transitionMatrix[]<-resEM$transitionMatrix
      model$emissionMatrix<-lapply(seq(dim(resEM$emissionArray)[3]), 
                                   function(i) resEM$emissionArray[ , 1:model$numberOfSymbols[i], i])
      names(model$emissionMatrix)<-model$channelNames
    }
    
  } else resEM <-NULL
  
  if(use.optimx){
    maxIP<-which.max(model$initialProbs)
    maxIPvalue<-model$initialProbs[maxIP]
    paramIP<-{function(x) which(x>0 & x<max(x))}(model$initialProbs)
    npIP<-length(paramIP)
    initNZ<-model$initialProbs>0
    initNZ[maxIP]<-2
    
    x<-which(model$transitionMatrix>0,arr.ind=TRUE)  
    transNZ<-x[order(x[,1]),]
    maxTM<-cbind(1:model$numberOfStates,max.col(model$transitionMatrix,ties.method="first")) #-1)*model$numberOfStates+1:model$numberOfStates
    maxTMvalue<-apply(model$transitionMatrix,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),]
    npTM<-nrow(paramTM)
    transNZ<-model$transitionMatrix>0
    transNZ[maxTM]<-2    
    
    
    if(model$numberOfChannels==1){
      
      x<-which(model$emissionMatrix>0,arr.ind=TRUE) 
      emissNZ<-x[order(x[,1]),]
      maxEM<-cbind(1:model$numberOfStates,max.col(model$emissionMatrix,ties.method="first"))
      maxEMvalue<-apply(model$emissionMatrix,1,max)
      paramEM <- rbind(emissNZ,maxEM)
      paramEM <- paramEM[!(duplicated(paramEM)|duplicated(paramEM,fromLast=TRUE)),]
      npEM<-nrow(paramEM)
      emissNZ<-model$emissionMatrix>0
      emissNZ[maxEM]<-2
      
      initialvalues<-c(
        if(npTM>0) model$transitionMatrix[paramTM],
        if(npEM>0) model$emissionMatrix[paramEM],
        if(npIP>0) model$initialProbs[paramIP]
      )
      initialvalues <- log(initialvalues)
      
      
      
      likfn<-function(pars,model,estimate=TRUE){
        pars<-exp(pars)         
        if(npTM>0){
          model$transitionMatrix[maxTM]<-maxTMvalue     
          model$transitionMatrix[paramTM]<-pars[1:npTM]
          model$transitionMatrix<-model$transitionMatrix/rowSums(model$transitionMatrix)  
        }
        if(npEM>0){
          model$emissionMatrix[maxEM]<-maxEMvalue     
          model$emissionMatrix[paramEM]<-pars[(npTM+1):(npTM+npEM)]
          model$emissionMatrix<-model$emissionMatrix/rowSums(model$emissionMatrix) 
        }
        
        if(npIP>0){
          model$initialProbs[maxIP]<-maxIPvalue
          model$initialProbs[paramIP]<-pars[(npTM+npEM+1):length(pars)]
          model$initialProbs[]<-model$initialProbs/sum(model$initialProbs)
        }
        
        if(estimate){
          - logLikHMM(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, obsArray)        
        } else model
      }
      
      rowSumsA<-rowSumsB<-rep(1,model$numberOfStates)
      sumInit<-1
      
      gradfn<-function(pars,model){
        pars<-exp(pars)         
        if(npTM>0){
          model$transitionMatrix[maxTM]<-maxTMvalue     
          model$transitionMatrix[paramTM]<-pars[1:npTM]
          rowSumsA<-rowSums(model$transitionMatrix)
          model$transitionMatrix<-model$transitionMatrix/rowSumsA          
        }
        if(npEM>0){
          model$emissionMatrix[maxEM]<-maxEMvalue     
          model$emissionMatrix[paramEM]<-pars[(npTM+1):(npTM+npEM)]
          rowSumsB<-rowSums(model$emissionMatrix)
          model$emissionMatrix<-model$emissionMatrix/rowSumsB          
        }
        
        if(npIP>0){
          model$initialProbs[maxIP]<-maxIPvalue
          model$initialProbs[paramIP]<-pars[(npTM+npEM+1):length(pars)]
          sumInit<-sum(model$initialProbs)
          model$initialProbs[]<-model$initialProbs/sumInit
        }         
        - gradient(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, obsArray, 
                   rowSumsA,rowSumsB,sumInit,transNZ,emissNZ,initNZ,pars) 
        
      }
    } else {      
      
      emissNZ<-lapply(model$emissionMatrix,function(i){
        x<-which(i>0,arr.ind=TRUE) 
        x[order(x[,1]),]
      })
      
      maxEM<-lapply(model$emissionMatrix,function(i) cbind(1:model$numberOfStates,max.col(i,ties.method="first")))
      
      maxEMvalue<-lapply(1:model$numberOfChannels, function(i) 
        apply(model$emissionMatrix[[i]],1,max))
      
      paramEM<-lapply(1:model$numberOfChannels,function(i) {
        x<-rbind(emissNZ[[i]],maxEM[[i]])
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),]
      })
      npEM<-sapply(paramEM,nrow)
      
      #       emissNZ<-sapply(1:model$numberOfChannels,function(i){
      #         x<-model$emissionMatrix[[i]]>0
      #         x[maxEM[[i]]]<-2
      #         x
      #       },simplify="array")  
      
      emissNZ<-array(0,c(model$numberOfStates,max(model$numberOfSymbols),model$numberOfChannels))
      for(i in 1:model$numberOfChannels){
        emissNZ[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]] > 0
        emissNZ[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-2
        
      }       
      
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
          model$transitionMatrix[maxTM]<-maxTMvalue     
          model$transitionMatrix[paramTM]<-pars[1:npTM]
          model$transitionMatrix<-model$transitionMatrix/rowSums(model$transitionMatrix)       
        }
        if(sum(npEM)>0){            
          for(i in 1:model$numberOfChannels){
            emissionArray[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]    
            emissionArray[,1:model$numberOfSymbols[i],i][paramEM[[i]]]<-pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])]
            rowSumsB<-rowSums(emissionArray[,1:model$numberOfSymbols[i],i])
            emissionArray[,1:model$numberOfSymbols[i],i]<-
              emissionArray[,1:model$numberOfSymbols[i],i]/rowSumsB
          }
        }
        
        if(npIP>0){
          model$initialProbs[maxIP]<-maxIPvalue
          model$initialProbs[paramIP]<-pars[(npTM+sum(npEM)+1):length(pars)]          
          model$initialProbs[]<-model$initialProbs/sum(model$initialProbs)
        }
        
        if(estimate){
          - logLikMCHMM(model$transitionMatrix, emissionArray, model$initialProbs, obsArray)
          
        } else {
          if(sum(npEM)>0){
            for(i in 1:model$numberOfChannels){
              model$emissionMatrix[[i]][]<-emissionArray[,1:model$numberOfSymbols[i],i]
            }
          }
          model
        }        
      }  
      rowSumsA<-rep(1,model$numberOfStates)
      rowSumsB<-matrix(1,model$numberOfStates,model$numberOfChannels)
      sumInit<-1
      
      gradfn<-function(pars,model){
        pars<-exp(pars)
        if(any(!is.finite(pars)) && estimate)
          return(.Machine$double.xmax)
        if(npTM>0){
          model$transitionMatrix[maxTM]<-maxTMvalue     
          model$transitionMatrix[paramTM]<-pars[1:npTM]
          rowSumsA<-rowSums(model$transitionMatrix)
          model$transitionMatrix<-model$transitionMatrix/rowSumsA         
        }
        if(sum(npEM)>0){            
          for(i in 1:model$numberOfChannels){
            emissionArray[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]    
            emissionArray[,1:model$numberOfSymbols[i],i][paramEM[[i]]]<-pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])]            
            rowSumsB[,i]<-rowSums(emissionArray[,1:model$numberOfSymbols[i],i])
            emissionArray[,1:model$numberOfSymbols[i],i]<-
              emissionArray[,1:model$numberOfSymbols[i],i]/rowSumsB[,i]
          }
        }
        
        if(npIP>0){
          model$initialProbs[maxIP]<-maxIPvalue
          model$initialProbs[paramIP]<-pars[(npTM+sum(npEM)+1):length(pars)]
          sumInit<-sum(model$initialProbs)
          model$initialProbs[]<-model$initialProbs/sumInit
        }        
        
        - gradientMC(model$transitionMatrix, emissionArray, model$initialProbs, obsArray,
                     rowSumsA,rowSumsB,sumInit,transNZ,emissNZ,initNZ,pars)
        
        
      }
    }
    
    
    if(is.null(optimx.control$fnscale) && use.em){
      optimx.control$fnscale <- -resEM$logLik 
    }
    if(is.null(optimx.control$kkt)){
      optimx.control$kkt<-FALSE
    }
    if(is.null(optimx.control$starttests)){
      optimx.control$starttests<-FALSE
    }
    resoptimx<-optimx(par=initialvalues,fn=likfn,gr=gradfn,method=method,itnmax=itnmax,control=optimx.control,model=model,...)
    model<-likfn(as.numeric(resoptimx[1:length(initialvalues)]),model,FALSE)
  } else resoptimx<-NULL
  
  list(model=model,logLik=ifelse(use.optimx,-resoptimx$value,resEM$logLik),em.result=resEM,optimx.result=resoptimx)
}