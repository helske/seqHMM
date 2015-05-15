#' Trim Small Probabilities of Hidden Markov Model
#' 
#' Function trimHMM tries to set small insignificant probabilities to zero 
#' without decreasing the likelihood significantly.
#' 
#' @export
#' @param model Model for which trimming is performed.
#' @param maxit Number of iterations. After zeroing small values, model is 
#'   refitted, and this is repeated until there is nothing to trim or maxit 
#'   iterations is used.
#' @param return.loglik Return the log-likelihood of trimmed model together with
#'   the model object. Default is FALSE.
#' @param zerotol Values smaller than this are trimmed to zero.
#' @param convergence.check Checks that the (local) maximum was reached when using \code{\link{optimx}} for model fitting. Defaults to \code{FALSE}. If model does not converge 
#' @param ... Further parameters passed on to \code{fitHMM}.
#'   
#' @seealso \code{\link{buildHMM}} for building Hidden Markov models before 
#'   fitting and \code{\link{fitHMM}} for fitting Hidden Markov models.
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
#' # Testing if changing parameter values smaller than 1e-03 to zero 
#' # leads to improved log-likelihood.
#' HMMtrim <- trimHMM(HMM$model, zerotol=1e-03, maxit=10)
#' 
trimHMM<-function(model,maxit=0,return.loglik=FALSE,zerotol=1e-8, 
                  convergence.check=FALSE,...){
  
  if(model$numberOfChannels==1){
    if(!(any(model$initialProbs<zerotol & model$initialProbs>0) || 
           any(model$transitionMatrix<zerotol & model$transitionMatrix>0)
         || any(model$emissionMatrix<zerotol & model$emissionMatrix>0))){

      if(convergence.check==TRUE){
        fit<-fitHMM(model, optimx.control=list(kkt=TRUE),...)
        if(fit$optimx.result$convcode==0 && fit$optimx.result$kkt1==TRUE && fit$optimx.result$kkt2==TRUE){
          print("Convergence check: (Local) optimum was found.")
        }else{
          print("Convergence check: Possible problem(s) with convergence.")
          if(fit$optimx.result$convcode!=0){
            print(paste("convcode =", fit$optimx.result$convcode))
          }
          if(fit$optimx.result$kkt1!=TRUE){
            print(paste("kkt1 =", fit$optimx.result$kkt1))
          }
          if(fit$optimx.result$kkt2!=TRUE){
            print(paste("kkt2 =", fit$optimx.result$kkt2))
          }
          print("Type help(optimx) for more information.")
        }
      }
      
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
    model$emissionMatrix[model$emissionMatrix<zerotol]<-0
    model$emissionMatrix<-model$emissionMatrix/rowSums(model$emissionMatrix)
    
    
    if(!is.finite(ll0<-logLik(model)))
      stop("Initial trimming resulted a non-finite log-likelihood. Try changing zerotol parameter.")
    
    if(maxit>0){
      for(ii in 1:maxit){
        fit<-fitHMM(model,...)
        ll<- fit$logLik
        
        if(ll>ll0){
          model<-fit$model
          ll0<-ll
        } else break
        
        if(!(any(model$initialProbs<zerotol & model$initialProbs>0) || 
               any(model$transitionMatrix<zerotol & model$transitionMatrix>0)
             || any(model$emissionMatrix<zerotol & model$emissionMatrix>0)))
          break
        
        model$initialProbs[model$initialProbs<zerotol]<-0
        model$initialProbs<-model$initialProbs/sum(model$initialProbs)
        model$transitionMatrix[model$transitionMatrix<zerotol]<-0
        model$transitionMatrix<-model$transitionMatrix/rowSums(model$transitionMatrix)
        model$emissionMatrix[model$emissionMatrix<zerotol]<-0
        model$emissionMatrix<-model$emissionMatrix/rowSums(model$emissionMatrix)
        
      }
      if(convergence.check==TRUE){
        fit<-fitHMM(model, optimx.control=list(kkt=TRUE),...)
        if(fit$optimx.result$convcode==0 && fit$optimx.result$kkt1==TRUE && fit$optimx.result$kkt2==TRUE){
          print("Convergence check: (Local) maximum was found.")
        }else{
          print("Convergence check: Possible problem(s) with convergence.")
          if(fit$optimx.result$convcode!=0){
            print(paste("convcode =", fit$optimx.result$convcode))
          }
          if(fit$optimx.result$kkt1!=TRUE){
            print(paste("kkt1 =", fit$optimx.result$kkt1))
          }
          if(fit$optimx.result$kkt2!=TRUE){
            print(paste("kkt2 =", fit$optimx.result$kkt2))
          }
          print("Type help(optimx) for more information.")
        }
      }
    }
    
  } else {
    if(!(any(model$initialProbs<zerotol & model$initialProbs>0) || 
           any(model$transitionMatrix<zerotol & model$transitionMatrix>0)
         || any(sapply(model$emissionMatrix,function(x) any(x<zerotol & x>0))))){
      if(convergence.check==TRUE){
        fit<-fitHMM(model, optimx.control=list(kkt=TRUE),...)
        if(fit$optimx.result$convcode==0 && fit$optimx.result$kkt1==TRUE && fit$optimx.result$kkt2==TRUE){
          print("Convergence check: (Local) optimum was found.")
        }else{
          print("Convergence check: Possible problem(s) with convergence.")
          if(fit$optimx.result$convcode!=0){
            print(paste("convcode =", fit$optimx.result$convcode))
          }
          if(fit$optimx.result$kkt1!=TRUE){
            print(paste("kkt1 =", fit$optimx.result$kkt1))
          }
          if(fit$optimx.result$kkt2!=TRUE){
            print(paste("kkt2 =", fit$optimx.result$kkt2))
          }
          print("Type help(optimx) for more information.")
        }
      }
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
    for(i in 1:model$numberOfChannels){
      model$emissionMatrix[[i]][model$emissionMatrix[[i]]<zerotol]<-0
      model$emissionMatrix[[i]]<-model$emissionMatrix[[i]]/
        rowSums(model$emissionMatrix[[i]])
      
    }
    
    if(!is.finite(ll0<-logLik(model)))
      stop("Initial trimming resulted a non-finite log-likelihood. Try changing zerotol parameter.")
    
    if(maxit>0){
      for(ii in 1:maxit){
        fit<-fitHMM(model,...)
        ll<- fit$logLik
        
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
        for(i in 1:model$numberOfChannels){
          model$emissionMatrix[[i]][model$emissionMatrix[[i]]<zerotol]<-0
          model$emissionMatrix[[i]]<-model$emissionMatrix[[i]]/
            rowSums(model$emissionMatrix[[i]])          
        }       
      }
    }
    if(convergence.check==TRUE){
      fit<-fitHMM(model, optimx.control=list(kkt=TRUE),...)
      if(fit$optimx.result$convcode==0 && fit$optimx.result$kkt1==TRUE && fit$optimx.result$kkt2==TRUE){
        print("Convergence check: (Local) optimum was found.")
      }else{
        print("Convergence check: Possible problem(s) with convergence.")
        if(fit$optimx.result$convcode!=0){
          print(paste("convcode =", fit$optimx.result$convcode))
        }
        if(fit$optimx.result$kkt1!=TRUE){
          print(paste("kkt1 =", fit$optimx.result$kkt1))
        }
        if(fit$optimx.result$kkt2!=TRUE){
          print(paste("kkt2 =", fit$optimx.result$kkt2))
        }
        print("Type help(optimx) for more information.")
      }
    }
  }
  if(maxit>0)
    print(paste(ii,"iteration(s) used."))
  
  if(ll0<ll_original){
    print(paste("Log-likelihood of the trimmed model is smaller than the original log-likelihood, ll_trim-ll_orig =", signif(ll0-ll_original, 3)))
  } else print(paste("Trimming improved log-likelihood, ll_trim-ll_orig =", signif(ll0-ll_original, 3)))
  
  
  
  if(return.loglik){
    list(model=model,loglik=ll0)
  } else model 
}