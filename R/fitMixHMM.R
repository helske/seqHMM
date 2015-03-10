#' Estimate Parameters of Mixture Hidden Markov Model

#' Function \code{fitMixHMM} estimates a mixture of hidden Markov models
#' using numerical maximization of log-likelihood. Initial values for estimation
#' are taken from the corresponding components of the model with preservation of
#' original zero probabilities.

#' @export
#' @import optimx
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
  
  stop("Under construction")
  
#   if(model$numberOfChannels==1){
#     # Using integers (Substract 1: indexing in C++ starts from 0)
#     obsArray <- data.matrix(model$observations)-1 
#     # Missing values
#     obsArray[obsArray>model$numberOfSymbols] <- model$numberOfSymbols
#     
#   } else{
#     obsArray<-array(0, c(model$numberOfSequences, model$lengthOfSequences, 
#                          model$numberOfChannels))
#     for(i in 1:model$numberOfChannels){
#       obsArray[,,i]<-data.matrix(model$observations[[i]])-1
#       obsArray[,,i][obsArray[,,i]>model$numberOfSymbols[i]] <- model$numberOfSymbols[i]
#     }           
#   }
#   #storage.mode(obsArray)<-"integer"
#   
#   
#   
#   # Largest initial probability
#   maxIP <- maxIPvalue <- npIP <- numeric(model$numberOfModels)
#   paramIP <- vector("list",model$numberOfModels)
#   for(m in 1:model$numberOfModels){
#     # Index of largest initial probability
#     maxIP[m] <- which.max(model$initialProbs[[m]])
#     # Value of largest initial probability
#     maxIPvalue[m] <- model$initialProbs[[m]][maxIP[m]]
#     # Rest of non-zero probs
#     paramIP[[m]] <- setdiff(which(model$initialProbs[[m]]>0),maxIP[m]) 
#     npIP[m] <- length(paramIP[[m]])
#   }
#   
#   # Largest transition probabilities (for each row)
#   transNZ <- maxTM <- maxTMvalue <- paramTM <- vector("list",model$numberOfModels)
#   npTM <- numeric(model$numberOfModels)
#   for(m in 1:model$numberOfModels){
#     x <- which(model$transitionMatrix[[m]]>0,arr.ind=TRUE)  
#     transNZ[[m]] <- x[order(x[,1]),]
#     maxTM[[m]] <- cbind(1:model$numberOfStates[[m]],
#                    max.col(model$transitionMatrix[[m]],ties.method="first"))
#     maxTMvalue[[m]] <- apply(model$transitionMatrix[[m]],1,max)
#     paramTM[[m]] <- rbind(transNZ[[m]],maxTM[[m]])
#     paramTM[[m]] <- paramTM[[m]][!(duplicated(paramTM[[m]])|duplicated(paramTM[[m]],fromLast=TRUE)),,drop=FALSE]
#     npTM[m] <- nrow(paramTM[[m]])
#   }
#   
#   # Single channel model
#   if(model$numberOfChannels==1){
#     emissNZ <- maxEM <- maxEMvalue <- paramEM <- vector("list",model$numberOfModels)
#     npEM <- numeric(model$numberOfModels)
#     # Largest emission probabilities (for each row)
#     for(m in 1:model$numberOfModels){
#       x <- which(model[[m]]$emissionMatrix>0,arr.ind=TRUE) 
#       emissNZ[[m]] <- x[order(x[,1]),]
#       maxEM[[m]] <- cbind(1:model[[m]]$numberOfStates,max.col(model[[m]]$emissionMatrix,ties.method="first"))
#       maxEMvalue[[m]] <- apply(model[[m]]$emissionMatrix,1,max)
#       paramEM[[m]] <- rbind(emissNZ[[m]],maxEM[[m]])
#       paramEM[[m]] <- paramEM[[m]][!(duplicated(paramEM[[m]])|duplicated(paramEM[[m]],fromLast=TRUE)),]
#       npEM[m] <- nrow(paramEM[[m]])
#     }
#     
#     # Initial parameter values (if anything to estimate)
#     initialvalues <- vector("list",model$numberOfModels)
#     for(m in 1:model$numberOfModels){
#       initialvalues[[m]] <- c(log(c(
#         if(npTM[[m]]>0) model$transitionMatrix[[m]][paramTM[[m]]],
#         if(npEM[[m]]>0) model$emissionMatrix[[m]][paramEM[[m]]],
#         if(npIP[[m]]>0) model$initialProbs[[m]][paramIP[[m]]]))
#       )
#     }
#     
#     
#     # Function for minimizing log likelihood
#     likfn<-function(pars,model,estimate=TRUE){
#       
#       for(m in 1:model$numberOfModels){
#         if(npTM[[m]]>0){
#           model$transitionMatrix[[m]][maxTM[[m]]] <- maxTMvalue[[m]] # Not needed?
#           # Exponentiate (need to be positive)
#           model$transitionMatrix[[m]][paramTM[[m]]] <- exp(pars[1:npTM[[m]]])
#           # Sum to 1
#           model$transitionMatrix[[m]] <- model$transitionMatrix[[m]]/
#             rowSums(model$transitionMatrix[[m]])  
#         }
#         if(npEM[[m]]>0){
#           model$emissionMatrix[[m]][maxEM[[m]]] <- maxEMvalue[[m]] # Not needed? 
#           model$emissionMatrix[[m]][paramEM[[m]]] <- exp(pars[(npTM[[m]]+1):
#                                                               (npTM[[m]]+npEM[[m]])])
#           model$emissionMatrix[[m]] <- model$emissionMatrix[[m]]/
#             rowSums(model$emissionMatrix[[m]]) 
#         }
#         
#         if(npIP[[m]]>0){
#           model$initialProbs[[m]][maxIP[[m]]] <- maxIPvalue[[m]] # Not needed?
#           model$initialProbs[[m]][paramIP[[m]]] <- exp(pars[npTM[[m]]+npEM[[m]]+
#                                                             1:npIP[[m]]])
#           model$initialProbs[[m]][] <- model$initialProbs[[m]]/
#             sum(model$initialProbs[[m]])
#         } 
#         
# #         if()
#           
#           
#           if(estimate){
#             # C++ function
#             - logLikMixHMM(model$transitionMatrix, 
#                            cbind(model$emissionMatrix,1), # prob=1 for missing observations
#                            model$initialProbs, obsArray, lweights)     
#             # Return model parameters
#           } else model
#       }
#     }
#     
#     
#     # Same for multichannel model  
#   } else {      
#     
#     emissNZ <- maxEM <- maxEMvalue <- paramEM <- npEM <- vector("list",model$numberOfModels)
# 
#     for(m in 1:model$numberOfModels){
#       emissNZ[[m]] <- lapply(model$emissionMatrix[[m]],function(i){
#         x<-which(i>0,arr.ind=TRUE) 
#         x[order(x[,1]),]
#       })
#       
#       maxEM[[m]] <- lapply(model$emissionMatrix[[m]],function(i) cbind(1:model$numberOfStates[[m]],max.col(i,ties.method="first")))
#       
#       maxEMvalue[[m]] <- lapply(1:model$numberOfChannels, function(i) 
#         apply(model$emissionMatrix[[m]][[i]],1,max))
#       
#       paramEM[[m]] <- lapply(1:model$numberOfChannels,function(i) {
#         x<-rbind(emissNZ[[m]][[i]],maxEM[[m]][[i]])
#         x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),,drop = FALSE]
#       })
#       npEM[[m]] <- sapply(paramEM[[m]],nrow)    
#     }
#     
# 
#     initialvalues <- vector("list",model$numberOfModels)
#     for(m in 1:model$numberOfModels){
#       initialvalues[[m]] <- c(log(c(
#         if(npTM[[m]]>0) model$transitionMatrix[[m]][paramTM[[m]]],
#         if(sum(npEM[[m]])>0) unlist(sapply(1:model$numberOfChannels,
#                                       function(x) model$emissionMatrix[[m]][[x]][paramEM[[m]][[x]]])),
#         if(npIP[[m]]>0) model$initialProbs[[m]][paramIP[[m]]]))
#       )
#     }
#     npModels <- sapply(initialvalues, length)
#     initialvalues <- unlist(initialvalues)
#     
#     emissionArray <- vector("list",model$numberOfModels)
#     for(m in 1:model$numberOfModels){
#       emissionArray[[m]] <- array(1,c(model$numberOfStates[[m]],
#                                       max(model$numberOfSymbols)+1,
#                                       model$numberOfChannels))
#       for(i in 1:model$numberOfChannels)
#         emissionArray[[m]][,1:model$numberOfSymbols[i],i] <- model$emissionMatrix[[m]][[i]]          
#     }
#     
#     # Tähän asti ok(?)
#     
#     likfn<-function(pars,model,estimate=TRUE){
#       
#       if(any(!is.finite(pars)) && estimate)
#         return(.Machine$double.xmax)
#       
#       
#       for(m in 1:model$numberOfModels){
#         csnpM <- cumsum(npModels[1:(m-1)])
#         if(npTM[[m]]>0){
#           model$transitionMatrix[[m]][maxTM[[m]]] <- maxTMvalue[[m]]     
#           model$transitionMatrix[[m]][paramTM[[m]]] <- 
#             exp(pars[csnpM+1:npTM[[m]]])
#           model$transitionMatrix[[m]] <- model$transitionMatrix[[m]]/
#             rowSums(model$transitionMatrix[[m]])       
#         }
#         if(sum(npEM[[m]])>0){            
#           for(i in 1:model$numberOfChannels){
#             emissionArray[[m]][,1:model$numberOfSymbols[i],i][maxEM[[m]][[i]]] <- 
#               maxEMvalue[[m]][[i]]    
#             emissionArray[[m]][,1:model$numberOfSymbols[i],i][paramEM[[m]][[i]]] <-
#               exp(pars[(csnpM+npTM[[m]]+1+c(0,cumsum(npEM[[m]]))[i]):
#                          (csnpM+npTM[[m]]+cumsum(npEM[[m]])[i])])
#             rowSumsB <- rowSums(emissionArray[[m]][,1:model$numberOfSymbols[i],i])
#             emissionArray[[m]][,1:model$numberOfSymbols[i],i] <-
#               emissionArray[[m]][,1:model$numberOfSymbols[i],i]/rowSumsB
#           }
#         }
#         
#         if(npIP[[m]]>0){
#           model$initialProbs[[m]][maxIP[[m]]] <- maxIPvalue[[m]]
#           model$initialProbs[[m]][paramIP[[m]]] <- 
#             exp(pars[(csnpM+npTM[[m]]+sum(npEM[[m]])+1):
#                        (csnpM+npTM[[m]]+sum(npEM[[m]])+npIP[[m]])])     
#           model$initialProbs[[m]][] <- model$initialProbs[[m]]/
#             sum(model$initialProbs[[m]])
#         } 
#       }
#       if(estimate){
#         - logLikMCHMM(model$transitionMatrix, emissionArray, model$initialProbs, 
#                       obsArray)     
#       } else {
#         for(m in 1:model$numberOfModels){
#           if(sum(npEM[[m]])>0){
#             for(i in 1:model$numberOfChannels){
#               model$emissionMatrix[[m]][[i]][] <- 
#                 emissionArray[[m]][,1:model$numberOfSymbols[i],i]
#             }
#           }
#         }
#         model
#       }        
#     }  
#     
#   }
#   
#   if(is.null(optimx.control$kkt)){
#     optimx.control$kkt <- FALSE
#   }
#   if(is.null(optimx.control$starttests)){
#     optimx.control$starttests <- FALSE
#   }
#   
#   resoptimx <- optimx(par=initialvalues, fn=likfn, method=method, 
#                       itnmax=itnmax, control=optimx.control, model=model,...)
#   model <- likfn(as.numeric(resoptimx[1:length(initialvalues)]), model, FALSE)
#   
#   
#   
#   list(model=model,logLik=-resoptimx$value,optimx.result=resoptimx)
}