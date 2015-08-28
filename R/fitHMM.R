#' Estimate Parameters of Hidden Markov Model
#'
#' Function \code{fitHMM} estimates the initial state, transition and emission 
#' probabilities of hidden Markov model. Initial values for estimation are taken from the 
#' corresponding components of the model with preservation of original zero 
#' probabilities.
#' 
#' By default, estimation start with EM algorithm and then switches to direct 
#' numerical maximization.
#' 
#' @export 
#' @import nloptr
#' @param model Hidden Markov model of class \code{HMModel}.
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
#'   \item{print_level}{Level of printing. Possible values are 0 
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
#'    \item{population}{\code{4*length(initialvalues)}} (number of starting points) 
#'}
#' @param lb,ub Lower and upper bounds for parameters in Softmax parameterization. 
#' Default interval is [pmin(-10,2*initialvalues), pmax(10,2*initialvalues)]. 
#' Used only in the global optimization step.
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
#' @seealso \code{\link{buildHMM}} for building Hidden Markov models before 
#'   fitting, \code{\link{trimHMM}} for finding better models by changing small
#'   parameter values to zero, \code{\link{BIC.HMModel}} for computing the
#'   value of the Bayesian information criterion of the model, and 
#'   \code{\link{plot.HMModel}} and \code{\link{ssplot}} for plotting 
#'   HMModel objects.
#' @details By default the \code{fitHMM} function uses only the \code{nloptr} function which 
#'   uses the multilevel single linkage method for global optimization 
#'   (\code{NLOPT_GD_MLSL} as \code{algorithm} in \code{global_control}). It performs 
#'   a sequence of local optimizations from random starting points, by default using 
#'   the BFGS algorithm (\code{NLOPT_LD_LBFGS} as \code{local_opts} in 
#'   \code{global_control}). The user can set the maximum number of evaluations or 
#'   limit the time used for the optimization.
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Initial values for emission matrices
#' B_child <- matrix(NA, nrow = 3, ncol = 2)
#' B_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 11:15])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow = 3, ncol = 2)
#' B_marr[1,] <- seqstatf(marr.seq[, 1:5])[, 2] + 0.1
#' B_marr[2,] <- seqstatf(marr.seq[, 6:10])[, 2] + 0.1
#' B_marr[3,] <- seqstatf(marr.seq[, 11:15])[, 2] + 0.1
#' B_marr <- B_marr / rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow = 3, ncol = 2)
#' B_left[1,] <- seqstatf(left.seq[, 1:5])[, 2] + 0.1
#' B_left[2,] <- seqstatf(left.seq[, 6:10])[, 2] + 0.1
#' B_left[3,] <- seqstatf(left.seq[, 11:15])[, 2] + 0.1
#' B_left <- B_left / rowSums(B_left)
#' 
#' # Initial values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow=3, ncol=3, byrow=TRUE)
#' 
#' # Initial values for initial state probabilities
#' initialProbs <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- buildHMM(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transitionMatrix = A,
#'   emissionMatrix = list(B_child, B_marr, B_left), 
#'   initialProbs = initialProbs
#'   )
#' 
#' # Fitting the model with different settings
#' 
#' # Only EM with default values
#' HMM1 <- fitHMM(bHMM, em_step = TRUE, global_step = FALSE, local_step = FALSE)
#' HMM1$logLik #-5507.003
#' 
#' \dontrun{
#' 
#' # EM with LBFGS
#' HMM2 <- fitHMM(bHMM, em_step = TRUE, global_step = FALSE, local_step = TRUE)
#' HMM2$logLik # -5507.003
#' 
#' # Only LBFGS
#' HMM3 <- fitHMM(bHMM, em_step = FALSE, global_step = FALSE, local_step = TRUE)
#' HMM3$logLik #-5493.271
#' 
#' # Global optimization via MLSL_LDS with LBFGS as local optimizer and final polisher
#' HMM4 <- fitHMM(bHMM, em_step = FALSE, global_step = TRUE, local_step = TRUE, 
#'   global_control = list(maxeval = 3000, maxtime = 0))
#' HMM4$logLik #-5403.383
#' 
#' # As previously, but now we use five iterations from EM algorithm for defining initial values and boundaries
#' # Note smaller maxeval for global optimization
#' HMM5 <- fitHMM(bHMM, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
#'   em_control = list(maxeval = 5), global_control = list(maxeval = 750, maxtime = 0))
#' HMM5$logLik #-5403.383
#' 
#' }
#' 
fitHMM<-function(model, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
  em_control=list(), global_control=list(), 
  local_control=list(), lb, ub, soft = TRUE, ...){
  
  if(!em_step && !global_step && !local_step){
    stop("No method chosen for estimation. Choose at least one from em_step, global_step, and local_step.")
  }
  
  if(model$numberOfChannels == 1){
    model$observations <- list(model$observations)
    model$emissionMatrix <- list(model$emissionMatrix)
  }
  
  obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
  for(i in 1:model$numberOfChannels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$numberOfSymbols[i]]<-model$numberOfSymbols[i]
  }
  
  if(em_step){
    em.con <- list(print_level = 0, maxeval = 100, reltol = 1e-8)
    nmsC <- names(em.con)  
    em.con[(namc <- names(em_control))] <- em_control
    if (length(noNms <- namc[!namc %in% nmsC])) 
      warning("Unknown names in em_control: ", paste(noNms, collapse = ", "))
    
    
    emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    if(soft){
      resEM<-EM(model$transitionMatrix, emissionArray, model$initialProbs, obsArray, 
        model$numberOfSymbols, em.con$maxeval, em.con$reltol,em.con$print_level)
    } else {
      resEM<-hardEMMC(model$transitionMatrix, emissionArray, model$initialProbs, obsArray, 
        model$numberOfSymbols, em.con$maxeval, em.con$reltol,em.con$print_level)
    }
    if(resEM$change< -1e-5)
      warning("EM algorithm stopped due to the decreasing log-likelihood. ")
    
    
    for(i in 1:model$numberOfChannels)
      model$emissionMatrix[[i]][]<-resEM$emissionArray[ , 1:model$numberOfSymbols[i], i]                                     
    
    
    model$initialProbs[]<-resEM$initialProbs
    model$transitionMatrix[]<-resEM$transitionMatrix
    ll <- resEM$logLik
  } else resEM <-NULL
  
  if(global_step || local_step){
    maxIP<-which.max(model$initialProbs)
    maxIPvalue<-model$initialProbs[maxIP]
    paramIP<-setdiff(which(model$initialProbs>0),maxIP)
    
    npIP<-length(paramIP)
    initNZ<-model$initialProbs>0
    initNZ[maxIP]<-0
    
    
    x<-which(model$transitionMatrix>0,arr.ind=TRUE)  
    transNZ<-x[order(x[,1]),]
    maxTM<-cbind(1:model$numberOfStates,max.col(model$transitionMatrix,ties.method="first"))
    maxTMvalue<-apply(model$transitionMatrix,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),,drop=FALSE]
    npTM<-nrow(paramTM)
    transNZ<-model$transitionMatrix>0
    transNZ[maxTM]<-0  
    
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
      emissNZ[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-0
      
    }       
    
    initialvalues<-c(log(c(
      if(npTM>0) model$transitionMatrix[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$numberOfChannels,
        function(x) model$emissionMatrix[[x]][paramEM[[x]]])),
      if(npIP>0) model$initialProbs[paramIP]))
    )         
    
    emissionArray<-array(1,c(model$numberOfStates,max(model$numberOfSymbols)+1,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]          
    
    objectivef<-function(pars, model, estimate = TRUE){
      
      if(any(!is.finite(exp(pars))) && estimate)
        return(.Machine$double.xmax^075)
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
          emissionArray[,1:model$numberOfSymbols[i],i]<-
            emissionArray[,1:model$numberOfSymbols[i],i]/rowSums(emissionArray[,1:model$numberOfSymbols[i],i, drop = FALSE])
        }
      }
      
      if(npIP>0){
        model$initialProbs[maxIP]<-maxIPvalue
        model$initialProbs[paramIP]<-exp(pars[(npTM+sum(npEM)+1):(npTM+sum(npEM)+npIP)])
        model$initialProbs[]<-model$initialProbs/sum(model$initialProbs)
      } 
      if(estimate){
        objective(model$transitionMatrix, emissionArray, model$initialProbs, obsArray,
          transNZ, emissNZ, initNZ, model$numberOfSymbols)
      } else {
        if(sum(npEM)>0){
          for(i in 1:model$numberOfChannels){
            model$emissionMatrix[[i]][]<-emissionArray[,1:model$numberOfSymbols[i],i]
          }
        }
        model
      }
    }
    
    
    if(global_step){
      
      if(missing(lb)){
        lb <- -10
      }
      if(missing(ub)){
        ub <- 10
      }
      lb <- pmin(lb, 2*initialvalues)
      ub <- pmax(ub, 2*initialvalues)
      
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
      
      globalres <- nloptr(x0 = initialvalues, eval_f = objectivef, lb = lb, ub = ub,
        opts = global_control, model = model, estimate = TRUE, ...)
      initialvalues <- globalres$solution
      model <- objectivef(globalres$solution, model, FALSE)
      ll <- -globalres$objective
    } else globalres <- NULL
    
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
        eval_f = objectivef,
        opts = local_control, model = model, estimate = TRUE, ub = rep(300,length(initialvalues)), ...)
      model <- objectivef(localres$solution,model, FALSE)
      ll <- -localres$objective
    } else localres <- NULL
    
    
  } else globalres <- localres <- NULL
  
  if(model$numberOfChannels == 1){
    model$observations <- model$observations[[1]]
    model$emissionMatrix <- model$emissionMatrix[[1]]
  }
  list(model = model, logLik = ll, 
    em_results = resEM[4:6], global_results = globalres, local_results = localres)
}