#' Transform Multichannel Hidden Markov Model to Single Channel Representation
#' 
#' @export
#' @param model Object of class HMModel.
#' @param combine.missing Controls whether combined states are coded missing
#'   (NA) if some of the channels include missing information. Defaults to
#'   TRUE.
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
#' HMM <- fitHMM(bHMM, em.control=list(maxit=100,reltol=1e-8),
#' itnmax=10000, method="BFGS")
#' 
#' HMM_SC <- MCtoSC(HMM$model)
#' 
#' @seealso \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and 
#'   fitting Hidden Markov models.

MCtoSC<-function(model, combine.missing=TRUE){
  if(model$numberOfChannels==1)
    return(model)
  
  B<-matrix(0,model$numberOfStates,prod(model$numberOfSymbols))
  
  colnames(B)<-apply(
    expand.grid(lapply(model$emissionMatrix,colnames)),                
    1,paste0,collapse="/")
  rownames(B)<-rownames(model$emissionMatrix[[1]])
  for(i in 1:model$numberOfStates){
    B[i,]<-apply(expand.grid(lapply(model$emissionMatrix,function(x) x[i,])),1,prod)   
  }
  
  
  modelx<-model
  modelx$emissionMatrix<-B
  modelx$numberOfSymbols<-ncol(B)
  modelx$numberOfChannels<-as.integer(1)
  modelx$symbolNames<-colnames(B)
  modelx$channelNames<-NULL
  
  modelx$observations<-model$observations[[1]]
  for(i in 2:model$numberOfChannels)
    modelx$observations<-as.data.frame(mapply(paste, modelx$observations,
                                              model$observations[[i]],
                                              USE.NAMES=FALSE,SIMPLIFY=FALSE,
                                              MoreArgs=list(sep="/")))
  names(modelx$observations)<-names(model$observations[[1]])   
  if(combine.missing==TRUE){
    modelx$observations[Reduce("|",
                               lapply(
                                 model$observations, 
                                 function(x) 
                                   x==attr(model$observations[[1]], "nr") |
                                   x==attr(model$observations[[1]], "void") |
                                   is.na(x)))]<-NA
  }
  modelx$observations<-seqdef(modelx$observations)
  attr(modelx$observations, "xtstep") <- attr(model$observations[[1]], "xtstep")
  attr(modelx$observations, "missing.color") <- attr(model$observations[[1]], "missing.color")
  attr(modelx$observations, "nr") <- attr(model$observations[[1]], "nr")
  attr(modelx$observations, "void") <- attr(model$observations[[1]], "void")
  attr(modelx$observations, "missing") <- attr(model$observations[[1]], "missing")
  attr(modelx$observations, "start") <- attr(model$observations[[1]], "start")
  modelx
}


# Jos logLik(modelx)!=logLik(model), pit?isi johtua puuttuvista
# (toinen sallii osittain puuttuvan tiedon):
# modelz<-model
# modelz$observations[[1]][is.na(fit$model$obs[[1]]) | 
#                           is.na(fit$model$obs[[2]]) | 
#                           is.na(fit$model$obs[[3]])]<-NA
# modelz$observations[[2]][is.na(fit$model$obs[[1]]) | 
#                           is.na(fit$model$obs[[2]]) | 
#                           is.na(fit$model$obs[[3]])]<-NA
# modelz$observations[[3]][is.na(fit$model$obs[[1]]) | 
#                           is.na(fit$model$obs[[2]]) | 
#                           is.na(fit$model$obs[[3]])]<-NA
# logLik(modelx)-logLik(modelz)