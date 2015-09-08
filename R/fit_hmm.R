#' Estimate Parameters of Hidden Markov Model
#'
#' Function \code{fit_hmm} estimates the initial state, transition and emission 
#' probabilities of hidden Markov model. Initial values for estimation are taken from the 
#' corresponding components of the model with preservation of original zero 
#' probabilities. By default, the estimation start with EM algorithm and then switches 
#' to direct numerical maximization.
#' 
#' @export 
#' @import nloptr
#' @param model Hidden Markov model of class \code{hmm}.
#' @param em_step Logical, use EM algorithm at the start of parameter estimation.
#'   The default is \code{TRUE}. Note that EM algorithm is faster than direct numerical optimization, 
#'   but is even more prone to get stuck in a local optimum.
#' @param global_step Logical, use global optimization via 
#'   \code{\link{nloptr}} (possibly after the EM step). The default is \code{TRUE}.
#'@param local_step Logical, use local optimization via 
#'   \code{\link{nloptr}} (possibly after the EM and/or global steps). The default is \code{TRUE}.
#' @param control_em Optional list of control parameters for for EM algorithm. 
#'   Possible arguments are \describe{ 
#'   \item{maxeval}{Maximum number of iterations, default is 100.} 
#'   \item{print_level}{Level of printing. Possible values are 0 
#'   (prints nothing), 1 (prints information at start and end of algorithm), and
#'   2 (prints at every iteration).} 
#'   \item{reltol}{Relative tolerance for convergence defined as \eqn{(sumLogLikNew - sumLogLikOld)/(abs(sumLogLikOld)+0.1)}. 
#'   Default is 1e-8.} }
#' @param control_global Optional list of additional arguments for 
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
#' @param control_local Optional list of additional arguments for 
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
#' @seealso \code{\link{build_hmm}} for building Hidden Markov models before 
#'   fitting, \code{\link{trim_hmm}} for finding better models by changing small
#'   parameter values to zero, \code{\link{BIC.hmm}} for computing the
#'   value of the Bayesian information criterion of the model, and 
#'   \code{\link{plot.hmm}} and \code{\link{ssplot}} for plotting 
#'   hmm objects.
#' @details The fitting function provides three estimation steps: 1) EM algorithm, 
#'   2) global optimization, and 3) local optimization. The user can call for one method 
#'   or any combination of these steps, but should note that they are preformed in the 
#'   above-mentioned order. The results from a former step are used as starting values 
#'   in a latter.
#' 
#'   By default the \code{fit_hmm} function starts with the EM algorithm,
#'   uses the multilevel single-linkage method (MLSL) with the LDS modification 
#'   for global optimization (\code{NLOPT_GD_MLSL_LDS} as \code{algorithm} in 
#'   \code{control_global}), and finishes with LBFGS as the local optimizer. 
#'   The MLSL method draws random starting points and performs a local optimization 
#'   from each. The LDS modification uses low-discrepancy sequences instead of 
#'   pseudo-random numbers as starting points and should improve the convergence rate. 
#'   By default, \code{fit_hmm} uses the BFGS algorithm as the local optimizer in the 
#'   MLSL (\code{NLOPT_LD_LBFGS} as \code{local_opts} in \code{control_global}). 
#'   In order to reduce the computation time spent on non-global optima, the 
#'   convergence tolerance of the local optimizer is set relatively large. At step 3, 
#'   a local optimization (BFGS by default) is run with a lower tolerance to find the 
#'   optimum with high precision.
#'   
#'   Any method available in the \code{nloptr} function can be used for the global and 
#'   local steps.
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' # Building one channel per type of event left, children or married
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
#' # Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Starting values for emission matrices
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
#' # Starting values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Starting values for initial state probabilities
#' init <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = init
#'   )
#' 
#' # Fitting the model with different settings
#' 
#' # Only EM with default values
#' HMM1 <- fit_hmm(bHMM, em_step = TRUE, global_step = FALSE, local_step = FALSE)
#' HMM1$logLik #-5507.003
#' 
#' \dontrun{
#' 
#' # EM with LBFGS
#' HMM2 <- fit_hmm(bHMM, em_step = TRUE, global_step = FALSE, local_step = TRUE)
#' HMM2$logLik # -5507.003
#' 
#' # Only LBFGS
#' HMM3 <- fit_hmm(bHMM, em_step = FALSE, global_step = FALSE, local_step = TRUE)
#' HMM3$logLik #-5493.271
#' 
#' # Global optimization via MLSL_LDS with LBFGS as local optimizer and final polisher
#' HMM4 <- fit_hmm(bHMM, em_step = FALSE, global_step = TRUE, local_step = TRUE, 
#'   control_global = list(maxeval = 3000, maxtime = 0))
#' HMM4$logLik #-5403.383
#' 
#' # As previously, but now we use five iterations from EM algorithm for defining initial values and boundaries
#' # Note smaller maxeval for global optimization
#' HMM5 <- fit_hmm(bHMM, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
#'   control_em = list(maxeval = 5), control_global = list(maxeval = 750, maxtime = 0))
#' HMM5$logLik #-5403.383
#' 
#' }
#' 
fit_hmm<-function(model, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
  control_em=list(), control_global=list(), 
  control_local=list(), lb, ub, ...){
  
  if(!em_step && !global_step && !local_step){
    stop("No method chosen for estimation. Choose at least one from em_step, global_step, and local_step.")
  }
  
  if(model$n_channels == 1){
    model$observations <- list(model$observations)
    model$emission_matrix <- list(model$emission_matrix)
  }
  
  obsArray<-array(0,c(model$n_sequences,model$length_of_sequences,model$n_channels))
  for(i in 1:model$n_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$n_symbols[i]]<-model$n_symbols[i]
  }
  
  if(em_step){
    em.con <- list(print_level = 0, maxeval = 100, reltol = 1e-8)
    nmsC <- names(em.con)  
    em.con[(namc <- names(control_em))] <- control_em
    if (length(noNms <- namc[!namc %in% nmsC])) 
      warning("Unknown names in control_em: ", paste(noNms, collapse = ", "))
    
    
    emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
    for(i in 1:model$n_channels)
      emissionArray[,1:model$n_symbols[i],i]<-model$emission_matrix[[i]]
    
    resEM<-EM(model$transition_matrix, emissionArray, model$initial_probs, obsArray, 
        model$n_symbols, em.con$maxeval, em.con$reltol,em.con$print_level)

    if(resEM$change< -1e-5)
      warning("EM algorithm stopped due to the decreasing log-likelihood. ")
    
    
    for(i in 1:model$n_channels)
      model$emission_matrix[[i]][]<-resEM$emissionArray[ , 1:model$n_symbols[i], i]                                     
    
    
    model$initial_probs[]<-resEM$initialProbs
    model$transition_matrix[]<-resEM$transitionMatrix
    ll <- resEM$logLik
  } else resEM <-NULL
  
  if(global_step || local_step){
    maxIP<-which.max(model$initial_probs)
    maxIPvalue<-model$initial_probs[maxIP]
    paramIP<-setdiff(which(model$initial_probs>0),maxIP)
    
    npIP<-length(paramIP)
    initNZ<-model$initial_probs>0
    initNZ[maxIP]<-0
    
    
    x<-which(model$transition_matrix>0,arr.ind=TRUE)  
    transNZ<-x[order(x[,1]),]
    maxTM<-cbind(1:model$n_states,max.col(model$transition_matrix,ties.method="first"))
    maxTMvalue<-apply(model$transition_matrix,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),,drop=FALSE]
    npTM<-nrow(paramTM)
    transNZ<-model$transition_matrix>0
    transNZ[maxTM]<-0  
    
    emissNZ<-lapply(model$emission_matrix,function(i){
      x<-which(i>0,arr.ind=TRUE) 
      x[order(x[,1]),]
    })
    
    maxEM<-lapply(model$emission_matrix,function(i) cbind(1:model$n_states,max.col(i,ties.method="first")))
    
    maxEMvalue<-lapply(1:model$n_channels, function(i) 
      apply(model$emission_matrix[[i]],1,max))
    
    paramEM<-lapply(1:model$n_channels,function(i) {
      x<-rbind(emissNZ[[i]],maxEM[[i]])
      x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),,drop = FALSE]
    })
    npEM<-sapply(paramEM,nrow)
    
    emissNZ<-array(0,c(model$n_states,max(model$n_symbols),model$n_channels))
    for(i in 1:model$n_channels){
      emissNZ[,1:model$n_symbols[i],i]<-model$emission_matrix[[i]] > 0
      if(model$n_states == 1){
        emissNZ[,1:model$n_symbols[i],i][maxEM[[i]][2]]<-0
      } else emissNZ[,1:model$n_symbols[i],i][maxEM[[i]]]<-0      
    }       
    
    initialvalues<-c(if((npTM+sum(npEM)+npIPAll)>0) log(c(
      if(npTM>0) model$transition_matrix[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$n_channels,
        function(x) model$emission_matrix[[x]][paramEM[[x]]])),
      if(npIP>0) model$initial_probs[paramIP]))
    )         
    
    emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
    for(i in 1:model$n_channels)
      emissionArray[,1:model$n_symbols[i],i]<-model$emission_matrix[[i]]          
    
    objectivef<-function(pars, model, estimate = TRUE){
      
      if(any(!is.finite(exp(pars))) && estimate)
        return(.Machine$double.xmax^075)
      if(npTM>0){
        model$transition_matrix[maxTM]<-maxTMvalue     
        model$transition_matrix[paramTM]<-exp(pars[1:npTM])
        model$transition_matrix<-model$transition_matrix/rowSums(model$transition_matrix)         
      }
      if(sum(npEM)>0){            
        for(i in 1:model$n_channels){
          emissionArray[,1:model$n_symbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]    
          emissionArray[,1:model$n_symbols[i],i][paramEM[[i]]]<-
            exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])      
          emissionArray[,1:model$n_symbols[i],i]<-
            emissionArray[,1:model$n_symbols[i],i]/rowSums(emissionArray[,1:model$n_symbols[i],i, drop = FALSE])
        }
      }
      
      if(npIP>0){
        model$initial_probs[maxIP]<-maxIPvalue
        model$initial_probs[paramIP]<-exp(pars[(npTM+sum(npEM)+1):(npTM+sum(npEM)+npIP)])
        model$initial_probs[]<-model$initial_probs/sum(model$initial_probs)
      } 
      if(estimate){
        objective(model$transition_matrix, emissionArray, model$initial_probs, obsArray,
          transNZ, emissNZ, initNZ, model$n_symbols)
      } else {
        if(sum(npEM)>0){
          for(i in 1:model$n_channels){
            model$emission_matrix[[i]][]<-emissionArray[,1:model$n_symbols[i],i]
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
      
      if(is.null(control_global$maxeval)){
        control_global$maxeval <- 10000
      }
      if(is.null(control_global$maxtime)){
        control_global$maxtime <- 60
      }
      if(is.null(control_global$algorithm)){
        control_global$algorithm <- "NLOPT_GD_MLSL_LDS"
        control_global$local_opts <- list(algorithm = "NLOPT_LD_LBFGS",  xtol_rel = 1e-4)
        control_global$ranseed <- 123
        control_global$population <- 4*length(initialvalues)
      }
      
      globalres <- nloptr(x0 = initialvalues, eval_f = objectivef, lb = lb, ub = ub,
        opts = control_global, model = model, estimate = TRUE, ...)
      initialvalues <- globalres$solution
      model <- objectivef(globalres$solution, model, FALSE)
      ll <- -globalres$objective
    } else globalres <- NULL
    
    if(local_step){
      if(is.null(control_local$maxeval)){
        control_local$maxeval <- 10000
      }
      if(is.null(control_local$maxtime)){
        control_local$maxtime <- 60
      }
      if(is.null(control_local$algorithm)){
        control_local$algorithm <- "NLOPT_LD_LBFGS"
        control_local$xtol_rel <- 1e-8
      }
      ub <- rep(300,length(initialvalues))
      ub <- pmax(ub, 2*initialvalues)
      localres<-nloptr(x0 = initialvalues, 
        eval_f = objectivef,
        opts = control_local, model = model, estimate = TRUE, ub = ub, ...)
      model <- objectivef(localres$solution,model, FALSE)
      ll <- -localres$objective
    } else localres <- NULL
    
    
  } else globalres <- localres <- NULL
  
  if(model$n_channels == 1){
    model$observations <- model$observations[[1]]
    model$emission_matrix <- model$emission_matrix[[1]]
  }
  list(model = model, logLik = ll, 
    em_results = resEM[4:6], global_results = globalres, local_results = localres)
}