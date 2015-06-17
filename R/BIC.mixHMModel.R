#' Compute BIC for Mixture Hidden Markov Model
#' 
#' @import stats
#' @export
#' @rdname BIC
#' @param object Object of class \code{mixHMModel}.
#' @param trim Should trimming be performed before computing BIC? See
#'   \code{\link{trimHMM}} for details. Default is FALSE.
#' @param ... further parameters to \code{\link{trimHMM}}.
#'   
#' @seealso \code{\link{buildMixHMM}} and \code{\link{fitMixHMM}} for building and
#'   fitting mixture hidden Markov models, and \code{\link{trimHMM}} for finding better
#'   models by changing small parameter values to zero.
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
#' ## Initial values for emission matrices 
#' alphabet(child.seq) # Check for order of observed states
#' B1_child <- matrix(NA, nrow=4, ncol=2) 
#' B1_child[1,] <- c(10,1) # High prob. for childless, low for children
#' B1_child[2,] <- c(10,1)
#' B1_child[3,] <- c(10,1)
#' B1_child[4,] <- c(1,10) # Low prob. for childless, high for children
#' B1_child <- B1_child/rowSums(B1_child)
#' B1_child
#' 
#' alphabet(marr.seq)
#' B1_marr <- matrix(NA, nrow=4, ncol=2) 
#' B1_marr[1,] <- c(1,10)
#' B1_marr[2,] <- c(1,10)
#' B1_marr[3,] <- c(10,1)
#' B1_marr[4,] <- c(7,1)
#' B1_marr <- B1_marr/rowSums(B1_marr)
#' B1_marr
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(NA, nrow=4, ncol=2) 
#' B1_left[1,] <- c(1,10)
#' B1_left[2,] <- c(10,1)
#' B1_left[3,] <- c(10,1)
#' B1_left[4,] <- c(10,1)
#' B1_left <- B1_left/rowSums(B1_left)
#' B1_left
#' 
#' B2_child <- matrix(NA, nrow=4, ncol=2) 
#' B2_child[1,] <- c(10,1) 
#' B2_child[2,] <- c(10,1)
#' B2_child[3,] <- c(10,1) 
#' B2_child[4,] <- c(1,10)
#' B2_child <- B2_child/rowSums(B2_child)
#' B2_child
#' 
#' B2_marr <- matrix(NA, nrow=4, ncol=2) 
#' B2_marr[1,] <- c(1,10)
#' B2_marr[2,] <- c(5,1)
#' B2_marr[3,] <- c(10,1)
#' B2_marr[4,] <- c(7,1)
#' B2_marr <- B2_marr/rowSums(B2_marr)
#' B2_marr
#' 
#' B2_left <- matrix(NA, nrow=4, ncol=2) 
#' B2_left[1,] <- c(1,10)
#' B2_left[2,] <- c(1,10)
#' B2_left[3,] <- c(1,5)
#' B2_left[4,] <- c(1,5)
#' B2_left <- B2_left/rowSums(B2_left)
#' B2_left
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'              0,    0.9, 0.07, 0.03, 
#'              0,      0,  0.9,  0.1, 
#'                0,      0,    0,    1), 
#'              nrow=4, ncol=4, byrow=TRUE)
#' 
#' A2 <- matrix(c(0.94, 0.04, 0.01, 0.01,
#'                0,    0.94, 0.05, 0.01, 
#'                0,       0,  0.9,  0.1, 
#'                0,       0,    0,    1), 
#'              nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initialProbs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initialProbs2 <- c(0.95, 0.03, 0.01, 0.01)
#' 
#' # Birth cohort
#' cohort <- cut(biofam$birthyr, c(1912, 1935, 1945, 1957))
#' biofam$cohort <- factor(cohort, labels=c("1913-1935", "1936-1945", "1946-1957"))
#' 
#' # Setting initial values for parameters
#' bmHMM <- buildMixHMM(observations=list(child.seq, marr.seq, left.seq), 
#'                      transitionMatrix=list(A1,A2), 
#'                      emissionMatrix=list(list(B1_child, B1_marr, B1_left),
#'                                          list(B2_child, B2_marr, B2_left)),
#'                      initialProbs=list(initialProbs1, initialProbs2), 
#'                      formula=~sex*cohort, data=biofam)
#' 
#' # Fitting mixture hidden Markov model
#' HMM <- fitMixHMM(bmHMM)
#' 
#' # Computing the BIC of the model
#' BIC(HMM$model)
#' 
BIC.mixHMModel<-function(object,trim=FALSE,...){
  
  #if(expand)
  #  object<-MCtoSC(object)
  
  if(trim){
    trimmed<-trimHMM(object,return.loglik=TRUE,...)
    object<-trimmed$model
    loglik<-trimmed$loglik
  } else loglik<-logLik(object)
  
  
  -2*loglik + log(object$numberOfSequences*object$lengthOfSequences)*
    (sum(unlist(object$initialProbs)>0)+sum(unlist(object$transitionMatrix)>0)+
       sum(unlist(object$emissionMatrix)>0))
}