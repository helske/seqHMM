#' Most Probable Path of Hidden States of Hidden Markov Model given the
#' sequence.
#' 
#' Function \code{mostProbablePath} computes the most probable path of the
#' hidden states of the hidden Markov model given the observed sequence.
#' 
#' @export
#' @param model Hidden Markov model of class \code{HMModel} or \code{MCHMModel}.
#' 
#' @return List which contains the most probable path of states (mpp) given the
#'   observations and its log-probability (logP). In case of multiple
#'   observations, most probable path is computed independently for each
#'   sequence.
#'   
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
#' HMM <- fitHMM(bHMM, em.control=list(maxit=100,reltol=1e-8), itnmax=10000, method="BFGS")
#'   
#' # Computing the most probable paths 
#' mpp <- mostProbablePath(HMM$model)$mpp
#'   
#' @seealso \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and 
#'   fitting Hidden Markov models, and \code{\link{seqIplot}} and
#'   \code{\link{plot.MCSP}} for plotting most probable paths.
#'   

mostProbablePath<-function(model){

  model$initialProbs <- log(model$initialProbs)
  model$initialProbs[!is.finite(model$initialProbs)]<- -0.1*.Machine$double.xmax
  model$transitionMatrix <- log(model$transitionMatrix)
  model$transitionMatrix[!is.finite(model$transitionMatrix)]<- -0.1*.Machine$double.xmax
  if(model$numberOfChannels==1){
    
    model$emissionMatrix <- log(model$emissionMatrix)
    model$emissionMatrix[!is.finite(model$emissionMatrix)]<- -0.1*.Machine$double.xmax
    obsArray<-data.matrix(model$observations)-1
    obsArray[obsArray>model$numberOfSymbols]<-model$numberOfSymbols
    storage.mode(obsArray)<-"integer"
    if(inherits(model,"mixHMModel")){
      out<-viterbiMix(model$transitionMatrix, cbind(model$emissionMatrix,0), 
                   model$initialProbs, obsArray, model$beta, 
                   model$X, model$numberOfStates)
    } else{
    out<-viterbi(model$transitionMatrix, cbind(model$emissionMatrix,0), 
                 model$initialProbs, obsArray)
    }
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
    
    if(inherits(model,"mixHMModel")){
      out<-viterbiMixMC(model$transitionMatrix, emissionArray, 
                      model$initialProbs, obsArray, model$beta, 
                      model$X, model$numberOfStates)
    } else{
      out<-viterbiMC(model$transitionMatrix, emissionArray, 
                     model$initialProbs, obsArray)
    }

    
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
