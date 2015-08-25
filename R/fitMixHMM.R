#' Estimate Parameters of Mixture Hidden Markov Model
#' 
#' Function \code{fitMixHMM} estimates a mixture of hidden Markov models
#' using numerical maximization of log-likelihood. Initial values for estimation
#' are taken from the corresponding components of the model with preservation of
#' original zero probabilities.
#' 
#' @export
#' @importFrom Matrix .bdiag
#' @param model Hidden Markov model of class \code{mixHMModel}.
#' @param em_step Logical, use EM algorithm at the start of parameter estimation.
#'   The default is \code{TRUE}. Note that EM algorithm is faster than direct numerical optimization, 
#'   but is even more prone to get stuck in a local optimum.
#' @param global_step Logical, use global optimization via 
#'   \code{\link{nloptr}} (possibly after the EM step). The default is \code{TRUE}.
#'@param local_step Logical, use local optimization via 
#'   \code{\link{nloptr}} (possibly after the EM and/or global steps). The default is \code{TRUE}.
#' @param em_control Optional list of control parameters for for EM algorithm. 
#'   Possible arguments are \describe{ 
#'   \item{maxeval}{Maximum number of iterations, default is 100.} 
#'   \item{trace}{Level of printing. Possible values are 0 
#'   (prints nothing), 1 (prints information at start and end of algorithm), and
#'   2 (prints at every iteration).} 
#'   \item{reltol}{Relative tolerance for convergence defined as \eqn{(sumLogLikNew - sumLogLikOld)/(abs(sumLogLikOld)+0.1)}. 
#'   Default is 1e-8.} }
#' @param global_control Optional list of additional arguments for 
#'   \code{\link{nloptr}} argument \code{opts}. The default values are
#'   \describe{
#'    \item{algorithm}{\code{"NLOPT_GD_MLSL_LDS"}}
#'    \item{local_opts}{\code{list(algorithm = "NLOPT_LD_LBFGS",  xtol_rel = 1e-4)}}
#'    \item{ranseed}{\code{123}}
#'    \item{maxeval}{\code{10000} (maximum number of iterations in global optimization algorithm)}
#'    \item{maxtime}{\code{60} (maximum run time in seconds)}
#'}
#' @param lb,ub Lower and upper bounds for parameters in Softmax parameterization. 
#' Default interval is [pmin(-10,2*initialvalues), pmax(10,2*initialvalues)]. Used only in global optimization step.
#' @param local_control Optional list of additional arguments for 
#'   \code{\link{nloptr}} argument \code{opts}. The default values are
#'   \describe{
#'    \item{algorithm}{\code{"NLOPT_LD_LBFGS"}}
#'    \item{xtol_rel}{\code{1e-8}}
#'    \item{maxeval}{\code{10000} (maximum number of iterations)}
#'    \item{maxtime}{\code{60} (maximum run time in seconds)}
#'   }
#' @param ... Additional arguments to nloptr
#' @return List with components \item{model}{Estimated model. } 
#'   \item{logLik}{Log-likelihood of the estimated model. } 
#'   \item{em_results}{Results after the EM step. } 
#'   \item{global_results}{Results after the global step. }
#'   \item{local_results}{Results after the local step. }
#' @seealso \code{\link{buildMixHMM}} for building mixture HMM's; 
#' \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and
#'   fitting hidden Markov models without covariates; \code{\link{plot.mixHMModel}}
#'   for plotting \code{mixHMModel} objects and \code{\link{mssplot}} for plotting
#'   stacked sequence plots of \code{mixHMModel} objects.
#' @examples 
#' \dontrun{
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[complete.cases(biofam[c(2:4)]),]
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf==4 | bf==5 | bf==6
#' married <- bf == 2 | bf== 3 | bf==6
#' left <- bf==1 | bf==3 | bf==5 | bf==6 | bf==7
#' 
#' children[children==TRUE] <- "Children"
#' children[children==FALSE] <- "Childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf==7)>0 & rowSums(bf==5)>0) | 
#'             (rowSums(bf==7)>0 & rowSums(bf==6)>0),]
#' children[rownames(bf) %in% rownames(div) & bf==7] <- "Children"
#' 
#' married[married==TRUE] <- "Married"
#' married[married==FALSE] <- "Single"
#' married[bf==7] <- "Divorced"
#' 
#' left[left==TRUE] <- "Left home"
#' left[left==FALSE] <- "With parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf==7)>0 & rowSums(bf==2)>0 & rowSums(bf==3)==0 &  
#'           rowSums(bf==5)==0 & rowSums(bf==6)==0) | 
#'            (rowSums(bf==7)>0 & rowSums(bf==4)>0 & rowSums(bf==3)==0 &  
#'           rowSums(bf==5)==0 & rowSums(bf==6)==0),]
#' left[rownames(bf) %in% rownames(wp) & bf==7] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children, start=15)
#' marr.seq <- seqdef(married, start=15)
#' left.seq <- seqdef(left, start=15)
#' 
#' ## Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' B1_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow=4, ncol=2, byrow=TRUE)
#' 
#' alphabet(marr.seq)                      
#' B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                     nrow=4, ncol=3, byrow=TRUE)                   
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01, # High probability for having left home
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow=4, ncol=2, byrow=TRUE)
#' 
#' # Cluster 2
#' B2_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.01, 0.99), nrow=4, ncol=2, byrow=TRUE)
#'                      
#' B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.29, 0.7, 0.01),
#'                    nrow=4, ncol=3, byrow=TRUE)                   
#' 
#' B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow=4, ncol=2, byrow=TRUE) 
#' 
#' # Cluster 3
#' B3_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.01, 0.99), nrow=6, ncol=2, byrow=TRUE)
#' 
#' B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                    nrow=6, ncol=3, byrow=TRUE)                   
#' 
#' B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.50, 0.50,
#'                     0.01, 0.99,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow=6, ncol=2, byrow=TRUE) 
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'                  0,    0.9, 0.07, 0.03, 
#'                  0,      0,  0.9,  0.1, 
#'                  0,      0,    0,    1), 
#'              nrow=4, ncol=4, byrow=TRUE)
#' 
#' A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
#'                  0,  0.7,  0.1,   0.1, 0.05, 0.05,
#'                  0,    0, 0.85,  0.01,  0.1, 0.04,
#'                  0,    0,    0,   0.9, 0.05, 0.05,
#'                  0,    0,    0,     0,  0.9,  0.1,
#'                  0,    0,    0,     0,    0,    1), 
#'              nrow=6, ncol=6, byrow=TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initialProbs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initialProbs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Creating covariate swiss
#' biofam$swiss <- biofam$nat_1_02=="Switzerland"
#' biofam$swiss[biofam$swiss==TRUE] <- "Swiss"
#' biofam$swiss[biofam$swiss==FALSE] <- "Other"
#' 
#' # Build mixture HMM
#' bmHMM <- buildMixHMM(observations=list(child.seq, marr.seq, left.seq), 
#'                        transitionMatrix=list(A1,A1,A2), 
#'                        emissionMatrix=list(list(B1_child, B1_marr, B1_left),
#'                                            list(B2_child, B2_marr, B2_left),
#'                                            list(B3_child, B3_marr, B3_left)),
#'                        initialProbs=list(initialProbs1, initialProbs1,
#'                                          initialProbs2), 
#'                        formula=~sex*birthyr+sex*swiss, data=biofam,
#'                        clusterNames=c("Cluster 1", "Cluster 2", "Cluster 3"),
#'                        channelNames=c("Parenthood", "Marriage", "Left home"),
#'                        )
#' 
#' mHMM <- fitMixHMM(bmHMM)
#' 
#' # Coefficients of covariates
#' mHMM$model$beta
#' 
#' # Probabilities of belonging to each model for the first six subjects
#' head(mHMM$model$clusterProb)
#' }


fitMixHMM <- function(model, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
  em_control = list(), global_control = list(), local_control = list(), lb, ub, soft = TRUE, ...){
  
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
  
  if(em_step){
    em.con <- list(trace = 0, maxit=100,reltol=1e-8)
    nmsC <- names(em.con)  
    em.con[(namc <- names(em.control))] <- em.control
    if (length(noNms <- namc[!namc %in% nmsC])) 
      warning("Unknown names in em.control: ", paste(noNms, collapse = ", "))
    
    if(model$numberOfChannels==1){
      
      if(soft){
        resEM <- EMx(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, 
          obsArray, model$numberOfSymbols,  model$beta, model$X, model$numberOfStatesInClusters, 
          em.con$maxit, em.con$reltol,em.con$trace)
      } else {
        resEM <- hardEMx(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, 
          obsArray, model$numberOfSymbols,  model$beta, model$X, model$numberOfStatesInClusters, 
          em.con$maxit, em.con$reltol,em.con$trace)
      }
      if(resEM$change< -1e-5)
        warning("EM algorithm stopped due to the decreasing log-likelihood. ")      
      
      model$emissionMatrix[]<-resEM$emissionMatrix[,1:model$numberOfSymbols]
      
    } else {    
      
      emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
      for(i in 1:model$numberOfChannels)
        emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
      
      resEM <- EMMCx(model$transitionMatrix, emissionArray, model$initialProbs, obsArray, 
        model$numberOfSymbols, model$beta, model$X, model$numberOfStatesInClusters, em.con$maxit, em.con$reltol,em.con$trace)
      if(!is.null(resEM$error))
        stop("Initial values for beta resulted non-finite cluster probabilities.")
      if(resEM$change< -1e-5)
        warning("EM algorithm stopped due to decreasing log-likelihood. ")
      
      
      for(i in 1:model$numberOfChannels)
        model$emissionMatrix[[i]][]<-resEM$emissionArray[ , 1:model$numberOfSymbols[i], i]                                     
    }
    
    if(global_step){
      k <- 0
      for(m in 1:model$numberOfClusters){
        original_model$initialProbs[[m]] <- unname(resEM$initialProbs[(k+1):(k+model$numberOfStatesInClusters[m])])
        k <- sum(model$numberOfStatesInClusters[1:m])
      }
    } else model$initialProbs[] <- resEM$initialProbs
    model$transitionMatrix[]<-resEM$transitionMatrix
    model$beta[]<-resEM$beta
    ll <- resEM$logLik
  } else resEM <-NULL
  
  if(global_step){
    maxIP <- maxIPvalue <- npIP <- numeric(original_model$numberOfClusters)  
    paramIP <-  initNZ <-vector("list",original_model$numberOfClusters)
    for(m in 1:original_model$numberOfClusters){
      # Index of largest initial probability
      maxIP[m] <- which.max(original_model$initialProbs[[m]])
      # Value of largest initial probability
      maxIPvalue[m] <- original_model$initialProbs[[m]][maxIP[m]]
      # Rest of non-zero probs
      paramIP[[m]] <- setdiff(which(original_model$initialProbs[[m]]>0),maxIP[m])
      npIP[m] <- length(paramIP[[m]])
      initNZ[[m]]<-original_model$initialProbs[[m]]>0
      initNZ[[m]][maxIP[m]]<-2
    }
    initNZ<-unlist(initNZ)
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
        if(npIPAll>0) unlist(sapply(1:original_model$numberOfClusters,function(m)
          if(npIP[m]>0) original_model$initialProbs[[m]][paramIP[[m]]]))
      )),
        model$beta[,-1]
      )
      
      # Function for minimizing log likelihood
      likfn<-function(pars,model,estimate=TRUE){
        
        if(any(!is.finite(exp(pars))) && estimate)
          return(.Machine$double.xmax)
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
        
        for(m in 1:original_model$numberOfClusters){
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
          - sum(logLikMixHMM(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, obsArray,
            model$beta, model$X, model$numberOfStatesInClusters))   
        } else model
      }
      sumInit<-rep(1,original_model$numberOfClusters)    
      rowSumsA<-rowSumsB<-rep(1,model$numberOfStates)
      
      gradfn<-function(pars,model, estimate){      
        
        if(any(!is.finite(exp(pars))))
          return(.Machine$double.xmax)
        
        if(npTM>0){
          model$transitionMatrix[maxTM]<-maxTMvalue      # Not needed?
          # Exponentiate (need to be positive)
          model$transitionMatrix[paramTM]<-exp(pars[1:npTM])
          # Sum to 1
          rowSumsA<-rowSums(model$transitionMatrix) 
          model$transitionMatrix<-model$transitionMatrix/rowSumsA
        }
        if(npEM>0){
          model$emissionMatrix[maxEM]<-maxEMvalue     
          model$emissionMatrix[paramEM]<-exp(pars[(npTM+1):(npTM+npEM)])
          rowSumsB <- rowSums(model$emissionMatrix) 
          model$emissionMatrix<-model$emissionMatrix/rowSumsB
        }
        
        for(m in 1:original_model$numberOfClusters){
          if(npIP[m]>0){
            original_model$initialProbs[[m]][maxIP[[m]]] <- maxIPvalue[[m]] # Not needed?
            original_model$initialProbs[[m]][paramIP[[m]]] <- exp(pars[npTM+npEM+c(0,cumsum(npIP))[m]+
                1:npIP[m]])
            sumInit[m]<-sum(original_model$initialProbs[[m]])
            original_model$initialProbs[[m]][] <- original_model$initialProbs[[m]]/
              sumInit[m]
          }
        }
        model$initialProbs <- unlist(original_model$initialProbs)
        model$beta[,-1] <- pars[npTM+npEM+npIPAll+1:npBeta]
        
        - gradientx(model$transitionMatrix, cbind(model$emissionMatrix,1), model$initialProbs, 
          obsArray, rowSumsA, rowSumsB, 
          sumInit, transNZ, emissNZ, initNZ, exp(pars[1:(npTM+sum(npEM)+npIPAll)]), 
          model$beta, model$X, model$numberOfStatesInClusters)
        
        
        
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
        if(npIPAll>0) unlist(sapply(1:original_model$numberOfClusters,function(m)
          if(npIP[m]>0) original_model$initialProbs[[m]][paramIP[[m]]]))
      )),
        model$beta[,-1]
      )         
      
      emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
      for(i in 1:model$numberOfChannels)
        emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]          
      
      
      
      likfn<-function(pars,model,estimate=TRUE){
        
        if(any(!is.finite(exp(pars))) && estimate)
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
        for(m in 1:original_model$numberOfClusters){
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
          - sum(logLikMixMCHMM(model$transitionMatrix, emissionArray, model$initialProbs, obsArray,
            model$beta, model$X, model$numberOfStatesInClusters)) 
        } else {
          if(sum(npEM)>0){
            for(i in 1:model$numberOfChannels){
              model$emissionMatrix[[i]][]<-emissionArray[,1:model$numberOfSymbols[i],i]
            }
          }
          model
        }        
      } 
      
      sumInit<-rep(1,original_model$numberOfClusters)    
      rowSumsA<-rep(1,model$numberOfStates)
      rowSumsB<-matrix(1,model$numberOfStates,model$numberOfChannels)
      
      gradfn<-function(pars,model, estimate){      
        
        if(any(!is.finite(exp(pars))))
          return(.Machine$double.xmax)
        if(npTM>0){
          model$transitionMatrix[maxTM]<-maxTMvalue     
          model$transitionMatrix[paramTM]<-exp(pars[1:npTM])
          rowSumsA <- rowSums(model$transitionMatrix) 
          model$transitionMatrix<-model$transitionMatrix/rowSumsA      
        }
        if(sum(npEM)>0){            
          for(i in 1:model$numberOfChannels){
            emissionArray[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]    
            emissionArray[,1:model$numberOfSymbols[i],i][paramEM[[i]]]<-
              exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])
            rowSumsB[,i]<-rowSums(emissionArray[,1:model$numberOfSymbols[i],i])
            emissionArray[,1:model$numberOfSymbols[i],i]<-
              emissionArray[,1:model$numberOfSymbols[i],i]/rowSumsB[,i]
          }
        }
        for(m in 1:original_model$numberOfClusters){
          if(npIP[m]>0){
            original_model$initialProbs[[m]][maxIP[[m]]] <- maxIPvalue[[m]] # Not needed?
            original_model$initialProbs[[m]][paramIP[[m]]] <- exp(pars[npTM+sum(npEM)+c(0,cumsum(npIP))[m]+
                1:npIP[m]])
            sumInit[m]<-sum(original_model$initialProbs[[m]])
            original_model$initialProbs[[m]][] <- original_model$initialProbs[[m]]/
              sumInit[m]
          }
        }
        model$initialProbs <- unlist(original_model$initialProbs)
        model$beta[,-1] <- pars[npTM+sum(npEM)+npIPAll+1:npBeta]
        
        - gradientMCx(model$transitionMatrix, emissionArray, model$initialProbs, obsArray, rowSumsA, rowSumsB, 
          sumInit, transNZ, emissNZ, initNZ, exp(pars[1:(npTM+sum(npEM)+npIPAll)]), 
          model$beta, model$X, model$numberOfStatesInClusters)
        
        
        
      }
    }
    
    if(missing(lb)){
      lb <- rep(-10, length(initialvalues))
    }
    if(missing(ub)){
      ub <-  rep(10, length(initialvalues))
    }
    lb <- pmin(lb, 2*initialvalues)
    ub <- pmax(ub, 2*initialvalues)
    
    if(global_step){
      if(is.null(global_control$maxeval)){
        global_control$maxeval <- 10000
      }
      if(is.null(global_control$maxtime)){
        global_control$maxtime <- 60
      }
      if(is.null(global_control$algorithm)){
        global_control$algorithm <- "NLOPT_GD_MLSL_LDS"
        global_control$local_opts <- list(algorithm = "NLOPT_LD_LBFGS",  xtol_rel = 1e-4)
        global_control$ranseed <- 123
        global_control$population <- 4*length(initialvalues)
      }
      
      globalres<-nloptr(x0 = initialvalues, eval_f = likfn, eval_grad_f = gradfn, lb = lb, ub = ub,
        opts=global_control, model=model, estimate = TRUE, ...)
      initialvalues <- globalres$solution
      model<-likfn(globalres$solution,model, FALSE)
      ll <- -globalres$objective
    } globalres <- NULL
    
    if(local_step){
      if(is.null(local_control$maxeval)){
        local_control$maxeval <- 10000
      }
      if(is.null(local_control$maxtime)){
        local_control$maxtime <- 60
      }
      if(is.null(local_control$algorithm)){
        local_control$algorithm <- "NLOPT_LD_LBFGS"
        local_control$xtol_rel <- 1e-8
      }
      localres<-nloptr(x0 = initialvalues, 
        eval_f = likfn, eval_grad_f = gradfn,
        opts = local_control, model = model, estimate = TRUE, ...)
      model<-likfn(localres$solution,model, FALSE)
      ll <- -localres$objective
    } else localres <- NULL
    
    rownames(model$beta) <- colnames(model$X)
    colnames(model$beta) <- model$clusterNames
    
  } else globalres <- localres <- NULL
  
  pr <- exp(model$X%*%model$beta)
  model$clusterProbabilities <- pr/rowSums(pr)
  list(model = spreadModels(model), logLik = ll, 
    em_results=resEM[4:6], global_results = globalres, local_results = localres)
  
}