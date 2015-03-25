#' Estimate Parameters of Mixture Hidden Markov Model
#' 
#' Function \code{fitMixHMM} estimates a mixture of hidden Markov models
#' using numerical maximization of log-likelihood. Initial values for estimation
#' are taken from the corresponding components of the model with preservation of
#' original zero probabilities.
#' 
#' @export
#' @importFrom Matrix .bdiag
#' @param model Hidden Markov model of class HMModel or MCHMModel.
#' @param method Optimization method used by \code{optimx}. Default is 
#'   \code{"BFGS"}. Note that \code{fitHMM} uses Softmax parameterization so 
#'   unconstrained optimization methods are used.
#' @param itnmax Maximum number of iterations use by \code{optimx}. Default is 
#'   10000.
#' @param optimx.control Optional list of additional arguments for 
#'   \code{\link{optimx}} argument \code{control}. Note that default values for 
#'   \code{starttests} and \code{kkt} are set to \code{FALSE}, which differs 
#'   from the default behaviour of \code{optimx}. If EM algorithm is used, 
#'   \code{fnscale} is also set to current optimum (unless modified by user).
#' @param ... Additional arguments to optimx.
#' @return List with components \item{model}{Estimated model. } 
#'   \item{logLik}{Log-likelihood of the estimated model. } 
#'   \item{optimx.results}{Results from direct numerical optimization via 
#'   \code{\link{optimx}}. }
#' @seealso \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and
#'   fitting hidden Markov models without covariates.
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf==4 | bf==5 | bf==6
#' married <- bf == 2 | bf== 3 | bf==6
#' left <- bf==1 | bf==3 | bf==5 | bf==6
#' 
#' children[children==TRUE] <- "Children"
#' children[children==FALSE] <- "Childless"
#' 
#' married[married==TRUE] <- "Married"
#' married[married==FALSE] <- "Single"
#' 
#' left[left==TRUE] <- "Left home"
#' left[left==FALSE] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Initial values for emission matrices
#' B_child <- matrix(NA, nrow=3, ncol=2)
#' B_child[1,] <- seqstatf(child.seq[,1:5])[,2]+0.1
#' B_child[2,] <- seqstatf(child.seq[,6:10])[,2]+0.1
#' B_child[3,] <- seqstatf(child.seq[,11:15])[,2]+0.1
#' B_child <- B_child/rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow=3, ncol=2)
#' B_marr[1,] <- seqstatf(marr.seq[,1:5])[,2]+0.1
#' B_marr[2,] <- seqstatf(marr.seq[,6:10])[,2]+0.1
#' B_marr[3,] <- seqstatf(marr.seq[,11:15])[,2]+0.1
#' B_marr <- B_marr/rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow=3, ncol=2)
#' B_left[1,] <- seqstatf(left.seq[,1:5])[,2]+0.1
#' B_left[2,] <- seqstatf(left.seq[,6:10])[,2]+0.1
#' B_left[3,] <- seqstatf(left.seq[,11:15])[,2]+0.1
#' B_left <- B_left/rowSums(B_left)
#' 
#' # Initial values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#' 0,    0.9,  0.1,
#' 0,      0,    1), 
#' nrow=3, ncol=3, byrow=TRUE)
#' 
#' # Initial values for initial state probabilities
#' initialProbs <- c(0.9,0.09,0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- buildHMM(observations=list(child.seq, marr.seq, left.seq), 
#' transitionMatrix=A,
#' emissionMatrix=list(B_child, B_marr, B_left), 
#' initialProbs=initialProbs)
#' 
#' # Fitting hidden Markov model
#' HMM <- fitHMM(bHMM, em.control=list(maxit=100,reltol=1e-8),
#' itnmax=10000, method="BFGS")
#' 


fitMixHMM<-function(model,method="BFGS",itnmax=10000,optimx.control=list(),...){
  
  original_model <- model
  model <- combineModels(model)
  
  if(model$numberOfChannels==1){
    # Using integers (Substract 1: indexing in C++ starts from 0)
    obsArray <- data.matrix(model$observations)-1 
    # Missing values
    obsArray[obsArray>model$numberOfSymbols] <- model$numberOfSymbols
    
  } else{
    obsArray<-array(0, c(model$numberOfSequences, model$lengthOfSequences, 
                         model$numberOfChannels))
    for(i in 1:model$numberOfChannels){
      obsArray[,,i]<-data.matrix(model$observations[[i]])-1
      obsArray[,,i][obsArray[,,i]>model$numberOfSymbols[i]] <- model$numberOfSymbols[i]
    }           
  }
  
  #   # Index of largest initial probability
  #   maxIP<-which.max(model$initialProbs)
  #   # Value of largest initial probability
  #   maxIPvalue<-model$initialProbs[maxIP]
  #   # Rest of non-zero probs
  #   paramIP<-setdiff(which(model$initialProbs>0),maxIP)
  #   npIP<-length(paramIP)
  #   initNZ<-model$initialProbs>0
  #   initNZ[maxIP]<-2
  maxIP <- maxIPvalue <- npIP <- numeric(original_model$numberOfModels)
  paramIP <- vector("list",original_model$numberOfModels)
  for(m in 1:original_model$numberOfModels){
    # Index of largest initial probability
    maxIP[m] <- which.max(original_model$initialProbs[[m]])
    # Value of largest initial probability
    maxIPvalue[m] <- original_model$initialProbs[[m]][maxIP[m]]
    # Rest of non-zero probs
    paramIP[[m]] <- setdiff(which(original_model$initialProbs[[m]]>0),maxIP[m])
    npIP[m] <- length(paramIP[[m]])
  }
  npIPAll <- sum(unlist(npIP))
  # Largest transition probabilities (for each row)
  x<-which(model$transitionMatrix>0,arr.ind=TRUE)  
  transNZ<-x[order(x[,1]),]
  maxTM<-cbind(1:model$numberOfStates,max.col(model$transitionMatrix,ties.method="first"))
  maxTMvalue<-apply(model$transitionMatrix,1,max)
  paramTM <- rbind(transNZ,maxTM)
  paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),,drop=FALSE]
  npTM<-nrow(paramTM)
  transNZ<-model$transitionMatrix>0
  transNZ[maxTM]<-2    
  
  npBeta<-length(model$beta[,-1])
  model$beta[,1] <- 0
  # Single channel model
  if(model$numberOfChannels==1){
    # Largest emission probabilities (for each row)
    x<-which(model$emissionMatrix>0,arr.ind=TRUE) 
    emissNZ<-x[order(x[,1]),]
    maxEM<-cbind(1:model$numberOfStates,max.col(model$emissionMatrix,ties.method="first"))
    maxEMvalue<-apply(model$emissionMatrix,1,max)
    paramEM <- rbind(emissNZ,maxEM)
    paramEM <- paramEM[!(duplicated(paramEM)|duplicated(paramEM,fromLast=TRUE)),]
    npEM<-nrow(paramEM)
    emissNZ<-model$emissionMatrix>0
    emissNZ[maxEM]<-2
    
    
    
    # Initial parameter values (if anything to estimate)
    initialvalues<-c(log(c(
      if(npTM>0) model$transitionMatrix[paramTM],
      if(npEM>0) model$emissionMatrix[paramEM],
      if(npIPAll>0) unlist(sapply(1:original_model$numberOfModels,function(m)
        if(npIP[m]>0) original_model$initialProbs[[m]][paramIP[[m]]]))
      )),
      model$beta[,-1]
    )
    
    # Function for minimizing log likelihood
    likfn<-function(pars,model,estimate=TRUE){
      
      if(npTM>0){
        model$transitionMatrix[maxTM]<-maxTMvalue      # Not needed?
        # Exponentiate (need to be positive)
        model$transitionMatrix[paramTM]<-exp(pars[1:npTM])
        # Sum to 1
        model$transitionMatrix<-model$transitionMatrix/rowSums(model$transitionMatrix)  
      }
      if(npEM>0){
        model$emissionMatrix[maxEM]<-maxEMvalue     
        model$emissionMatrix[paramEM]<-exp(pars[(npTM+1):(npTM+npEM)])
        model$emissionMatrix<-model$emissionMatrix/rowSums(model$emissionMatrix) 
      }
      
      for(m in 1:original_model$numberOfModels){
        if(npIP[m]>0){
          original_model$initialProbs[[m]][maxIP[[m]]] <- maxIPvalue[[m]] # Not needed?
          original_model$initialProbs[[m]][paramIP[[m]]] <- exp(pars[npTM+npEM+c(0,cumsum(npIP))[m]+
                                                                       1:npIP[m]])
          original_model$initialProbs[[m]][] <- original_model$initialProbs[[m]]/
            sum(original_model$initialProbs[[m]])
        }
      }
      model$initialProbs <- unlist(original_model$initialProbs)
      model$beta[,-1] <- pars[npTM+npEM+npIPAll+1:npBeta]
      if(estimate){
        - logLikMixHMM(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, obsArray,
                       model$beta, model$X, original_model$numberOfStates)      
      } else model
    }
    # Same for multichannel model  
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
      x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),,drop = FALSE]
    })
    npEM<-sapply(paramEM,nrow)
    
    emissNZ<-array(0,c(model$numberOfStates,max(model$numberOfSymbols),model$numberOfChannels))
    for(i in 1:model$numberOfChannels){
      emissNZ[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]] > 0
      emissNZ[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-2
      
    }       
    
    initialvalues<-c(log(c(
      if(npTM>0) model$transitionMatrix[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$numberOfChannels,
                                    function(x) model$emissionMatrix[[x]][paramEM[[x]]])),
      if(npIPAll>0) unlist(sapply(1:original_model$numberOfModels,function(m)
        if(npIP[m]>0) original_model$initialProbs[[m]][paramIP[[m]]]))
      )),
      model$beta[,-1]
    )         
    
    emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]          
    
    
    
    likfn<-function(pars,model,estimate=TRUE){
      
      if(any(!is.finite(pars)) && estimate)
        return(.Machine$double.xmax)
      if(npTM>0){
        model$transitionMatrix[maxTM]<-maxTMvalue     
        model$transitionMatrix[paramTM]<-exp(pars[1:npTM])
        model$transitionMatrix<-model$transitionMatrix/rowSums(model$transitionMatrix)       
      }
      if(sum(npEM)>0){            
        for(i in 1:model$numberOfChannels){
          emissionArray[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]    
          emissionArray[,1:model$numberOfSymbols[i],i][paramEM[[i]]]<-
            exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])
          rowSumsB<-rowSums(emissionArray[,1:model$numberOfSymbols[i],i])
          emissionArray[,1:model$numberOfSymbols[i],i]<-
            emissionArray[,1:model$numberOfSymbols[i],i]/rowSumsB
        }
      }
      for(m in 1:original_model$numberOfModels){
        if(npIP[m]>0){
          original_model$initialProbs[[m]][maxIP[[m]]] <- maxIPvalue[[m]] # Not needed?
          original_model$initialProbs[[m]][paramIP[[m]]] <- exp(pars[npTM+sum(npEM)+c(0,cumsum(npIP))[m]+
                                                                       1:npIP[m]])
          original_model$initialProbs[[m]][] <- original_model$initialProbs[[m]]/
            sum(original_model$initialProbs[[m]])
        }
      }
      model$initialProbs <- unlist(original_model$initialProbs)
      model$beta[,-1] <- pars[npTM+sum(npEM)+npIPAll+1:npBeta]

      if(estimate){   
        - logLikMixMCHMM(model$transitionMatrix, emissionArray, model$initialProbs, obsArray,
                         model$beta, model$X, original_model$numberOfStates)   
      } else {
        if(sum(npEM)>0){
          for(i in 1:model$numberOfChannels){
            model$emissionMatrix[[i]][]<-emissionArray[,1:model$numberOfSymbols[i],i]
          }
        }
        model
      }        
    } 
  }
  
  if(is.null(optimx.control$kkt)){
    optimx.control$kkt <- FALSE
  }
  if(is.null(optimx.control$starttests)){
    optimx.control$starttests <- FALSE
  }
  resoptimx <- optimx(par=initialvalues, fn=likfn, method=method, 
                      itnmax=itnmax, control=optimx.control, model=model,...)
  model <- likfn(as.numeric(resoptimx[1:length(initialvalues)]), model, FALSE)
  
  list(model=(model),logLik=-resoptimx$value,optimx.result=resoptimx)
}