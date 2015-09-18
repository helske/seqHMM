#' Estimate Parameters of Mixture Hidden Markov Model
#' 
#' Function \code{fit_mhmm} estimates a mixture of hidden Markov models
#' using numerical maximization of log-likelihood. Initial values for estimation
#' are taken from the corresponding components of the model with preservation of
#' original zero probabilities.
#' 
#' @export
#' @param model Hidden Markov model of class \code{mhmm}.
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
#'   (prints nothing), 1 (prints information at start and end of algorithm), 
#'   2 (prints at every iteration), and 3 (prints also inside the coefficient optimization).} 
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
#'}
#' @param lb,ub Lower and upper bounds for parameters in Softmax parameterization. 
#' Default interval is [pmin(-10,2*initialvalues), pmax(10,2*initialvalues)]. Lower bound is used only in global optimization.
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
#' @seealso \code{\link{build_mhmm}} for building MHMMs; 
#' \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and
#'   fitting hidden Markov models; \code{\link{plot.mhmm}}
#'   for plotting \code{mhmm} objects and \code{\link{mssplot}} for plotting
#'   stacked sequence plots of \code{mhmm} objects.
#' @examples 
#' \dontrun{
#' require(TraMineR)
#' 
#' # Single-channel
#' 
#' data(mvad)
#' 
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school", 
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education", 
#'                  "joblessness", "school", "training")
#' mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, 
#'                    labels = mvad.labels, xtstep = 6)
#'                    
#' # Starting values for emission matrices
#' emiss_1 <- matrix(
#' c(0.26, 0.39, 0.01, 0.06, 0.04, 0.24,
#'   0.58, 0.12, 0.09, 0.10, 0.01, 0.10,
#'   0.73, 0.02, 0.09, 0.13, 0.01, 0.02), 
#' nrow = 3, ncol = 6, byrow = TRUE)
#' 
#' emiss_2 <- matrix(
#' c(0.01, 0.02, 0.01, 0.01, 0.94, 0.01,
#'   0.05, 0.06, 0.15, 0.01, 0.72, 0.01,
#'   0.19, 0.13, 0.60, 0.01, 0.05, 0.02,
#'   0.32, 0.03, 0.60, 0.03, 0.01, 0.01), 
#' nrow = 4, ncol = 6, byrow = TRUE)
#' 
#' # Starting values for transition matrices
#' 
#' tr_1 <-  matrix(
#'   c(0.80, 0.10, 0.10,
#'     0.10, 0.80, 0.10,
#'     0.10, 0.10, 0.80), 
#'   nrow=3, ncol=3, byrow=TRUE)
#' 
#' tr_2 <-  matrix(
#'   c(0.80, 0.10, 0.05, 0.05,
#'     0.05, 0.80, 0.10, 0.05,
#'     0.05, 0.05, 0.80, 0.10,
#'     0.05, 0.05, 0.10, 0.80), 
#'   nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Starting values for initial state probabilities
#' init_1 <- c(0.4, 0.3, 0.3)
#' init_2 <- c(0.3, 0.3, 0.2, 0.2)
#' 
#' # Building a MHMM with the starting values
#' bmhmm_mvad <- build_mhmm(
#'   observations = mvad.seq, 
#'   transition_matrix = list(tr_1, tr_2), 
#'   emission_matrix = list(emiss_1, emiss_2), 
#'   initial_probs = list(init_1, init_2)
#' )                   
#' 
#' # Fitting the model with the EM algorithm
#' mhmm_mvad <- fit_mhmm(
#'   bmhmm_mvad, em_step = TRUE, global_step = FALSE, local_step = FALSE)
#' 
#' ##############################################################
#' 
#' # Multichannel
#' 
#' data(biofam3c)
#' 
#' # Building sequence objects
#' child.seq <- seqdef(biofam3c$children)
#' marr.seq <- seqdef(biofam3c$married)
#' left.seq <- seqdef(biofam3c$left)
#' 
#' ## Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' emiss_1_child <- matrix(
#'   c(0.99, 0.01, # High probability for childless
#'     0.99, 0.01,
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' alphabet(marr.seq)                      
#' emiss_1_marr <- matrix(
#'   c(0.01, 0.01, 0.98, # High probability for single
#'     0.01, 0.01, 0.98,
#'     0.01, 0.98, 0.01, # High probability for married
#'     0.98, 0.01, 0.01), # High probability for divorced
#'   nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' alphabet(left.seq)
#' emiss_1_left <- matrix(
#'   c(0.01, 0.99, # High probability for living with parents
#'     0.99, 0.01, # High probability for having left home
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' # Cluster 2
#' emiss_2_child <- matrix(
#'   c(0.99, 0.01, # High probability for childless
#'     0.99, 0.01,
#'     0.99, 0.01,
#'     0.01, 0.99), 
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'                      
#' emiss_2_marr <- matrix(
#'   c(0.01, 0.01, 0.98, # High probability for single
#'     0.01, 0.01, 0.98,
#'     0.01, 0.98, 0.01, # High probability for married
#'     0.29, 0.7, 0.01),
#'   nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' emiss_2_left <- matrix(
#'   c(0.01, 0.99, # High probability for living with parents
#'     0.99, 0.01,
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 4, ncol = 2, byrow = TRUE) 
#' 
#' # Cluster 3
#' emiss_3_child <- matrix(
#'   c(0.99, 0.01, # High probability for childless
#'     0.99, 0.01,
#'     0.01, 0.99,
#'     0.99, 0.01,
#'     0.01, 0.99,
#'     0.01, 0.99), 
#'   nrow = 6, ncol = 2, byrow = TRUE)
#' 
#' emiss_3_marr <- matrix(
#'   c(0.01, 0.01, 0.98, # High probability for single
#'     0.01, 0.01, 0.98,
#'     0.01, 0.01, 0.98,
#'     0.01, 0.98, 0.01,
#'     0.01, 0.98, 0.01, # High probability for married
#'     0.98, 0.01, 0.01), # High probability for divorced
#'   nrow = 6, ncol = 3, byrow = TRUE)                   
#' 
#' emiss_3_left <- matrix(
#'   c(0.01, 0.99, # High probability for living with parents
#'     0.99, 0.01,
#'     0.50, 0.50,
#'     0.01, 0.99,
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 6, ncol = 2, byrow = TRUE) 
#' 
#' # Initial values for transition matrices
#' trans_1 <- matrix(
#'   c(0.80,   0.16, 0.03, 0.01,
#'        0,   0.90, 0.07, 0.03, 
#'        0,      0, 0.90, 0.10, 
#'        0,      0,    0,    1), 
#'   nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' trans_2 <- matrix(
#'   c(0.80, 0.10, 0.05,  0.03, 0.01, 0.01,
#'        0, 0.70, 0.10,  0.10, 0.05, 0.05,
#'        0,    0, 0.85,  0.01, 0.10, 0.04,
#'        0,    0,    0,  0.90, 0.05, 0.05,
#'        0,    0,    0,     0, 0.90, 0.10,
#'        0,    0,    0,     0,    0,    1), 
#'   nrow = 6, ncol = 6, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Birth cohort
#' biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
#' biofam3c$covariates$cohort <- factor(
#'   biofam3c$covariates$cohort, labels=c("1909-1935", "1936-1945", "1946-1957"))
#' 
#' # Build mixture HMM
#' bmhmm <- build_mhmm(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_matrix = list(trans_1, trans_1, trans_2),
#'   emission_matrix = list(list(emiss_1_child, emiss_1_marr, emiss_1_left),
#'                         list(emiss_2_child, emiss_2_marr, emiss_2_left), 
#'                         list(emiss_3_child, emiss_3_marr, emiss_3_left)),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   formula = ~sex + cohort, data = biofam3c$covariates,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Parenthood", "Marriage", "Left home"))
#' 
#' # Fitting the model with different settings
#' 
#' # Only EM with default values
#' mhmm_1 <- fit_mhmm(
#'   bmhmm, em_step = TRUE, global_step = FALSE, local_step = FALSE)
#' mhmm_1$logLik # -3081.383
#' 
#' \dontrun{
#' # EM with LBFGS
#' mhmm_2 <- fit_mhmm(
#'   bmhmm, em_step = TRUE, global_step = FALSE, local_step = TRUE)
#' mhmm_2$logLik # -3081.383
#' 
#' # Only LBFGS
#' mhmm_3 <- fit_mhmm(
#'   bmhmm, em_step = FALSE, global_step = FALSE, local_step = TRUE,
#'   control_local = list(maxeval = 5000, maxtime = 0))
#' mhmm_3$logLik # -3087.499373
#' 
#' # Global optimization via MLSL_LDS with LBFGS as local optimizer and final polisher
#' mhmm_4 <- fit_mhmm(
#'   bmhmm, em_step = FALSE, global_step = TRUE, local_step = TRUE, 
#'   control_global = list(maxeval = 5000, maxtime = 0))
#' mhmm_4$logLik # -3150.796
#' 
#' # As previously, but now we use ten iterations from EM algorithm for 
#' # defining initial values and boundaries
#' # Note smaller maxeval for global optimization
#' mhmm_5 <- fit_mhmm(
#'   bmhmm, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
#'   control_em = list(maxeval = 10), control_global = list(maxeval = 1000, maxtime = 0),
#'   control_local = list(maxeval = 500, maxtime = 0))
#' mhmm_5$logLik #-3081.383
#' }
#' 
#' # Coefficients of covariates
#' coef(mhmm_1$model)
#' 
#' # Probabilities of belonging to each model for the first six subjects
#' head(mhmm_1$model$cluster_prob)
#' }


fit_mhmm <- function(model, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
  control_em = list(), control_global = list(), control_local = list(), lb, ub, ...){
  
  if(!inherits(model, "mhmm"))
    stop("Argument model must be an object of class 'mhmm'.")
  
  df <- attr(model, "df")
  nobs <- attr(model, "nobs")
  original_model <- model
  model <- combine_models(model)
  
  if(model$n_channels == 1){
    model$observations <- list(model$observations)
    model$emission_matrix <- list(model$emission_matrix)
  }
  
  obsArray<-array(0, c(model$n_sequences, model$length_of_sequences, 
    model$n_channels))
  for(i in 1:model$n_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$n_symbols[i]] <- model$n_symbols[i]
  } 
  emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
  for(i in 1:model$n_channels)
    emissionArray[,1:model$n_symbols[i],i]<-model$emission_matrix[[i]]
  
  if(em_step){
    em.con <- list(print_level = 0, maxeval=100, reltol=1e-8)
    nmsC <- names(em.con)  
    em.con[(namc <- names(control_em))] <- control_em
    if (length(noNms <- namc[!namc %in% nmsC])) 
      warning("Unknown names in control_em: ", paste(noNms, collapse = ", "))
    
    resEM <- EMx(model$transition_matrix, emissionArray, model$initial_probs, obsArray, 
      model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters, em.con$maxeval, em.con$reltol,em.con$print_level)
    if(!is.null(resEM$error))
      stop("Initial values for coefficients of covariates resulted non-finite cluster probabilities.")
    if(resEM$change< -1e-5)
      warning("EM algorithm stopped due to decreasing log-likelihood. ")
    
    emissionArray <- resEM$emissionArray
    for(i in 1:model$n_channels)
      model$emission_matrix[[i]][]<-emissionArray[ , 1:model$n_symbols[i], i]
    
    if(global_step || local_step){
      k <- 0
      for(m in 1:model$n_clusters){
        original_model$initial_probs[[m]] <- unname(resEM$initialProbs[(k+1):(k+model$n_states_in_clusters[m])])
        k <- sum(model$n_states_in_clusters[1:m])
      }
    }
    
    model$transition_matrix[]<-resEM$transitionMatrix
    model$coefficients[]<-resEM$coefficients
    ll <- resEM$logLik
  } else resEM <-NULL
  
  if(global_step || local_step){
    maxIP <- maxIPvalue <- npIP <- numeric(original_model$n_clusters)  
    paramIP <-  initNZ <-vector("list",original_model$n_clusters)
    for(m in 1:original_model$n_clusters){
      # Index of largest initial probability
      maxIP[m] <- which.max(original_model$initial_probs[[m]])
      # Value of largest initial probability
      maxIPvalue[m] <- original_model$initial_probs[[m]][maxIP[m]]
      # Rest of non-zero probs
      paramIP[[m]] <- setdiff(which(original_model$initial_probs[[m]]>0),maxIP[m])
      npIP[m] <- length(paramIP[[m]])
      initNZ[[m]]<-original_model$initial_probs[[m]]>0
      initNZ[[m]][maxIP[m]]<-0
    }
    initNZ<-unlist(initNZ)
    npIPAll <- sum(unlist(npIP))
    # Largest transition probabilities (for each row)
    x<-which(model$transition_matrix>0,arr.ind=TRUE)  
    transNZ<-x[order(x[,1]),]
    maxTM<-cbind(1:model$n_states,max.col(model$transition_matrix,ties.method="first"))
    maxTMvalue<-apply(model$transition_matrix,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),,drop=FALSE]
    npTM<-nrow(paramTM)
    transNZ<-model$transition_matrix>0
    transNZ[maxTM]<-0    
    
    npCoef<-length(model$coefficients[,-1])
    model$coefficients[,1] <- 0
    
    
    emissNZ<-lapply(model$emission_matrix,function(i){
      x<-which(i>0,arr.ind=TRUE) 
      x[order(x[,1]),]
    })
    
    if(model$n_states > 1){
      maxEM <- lapply(model$emission_matrix,function(i) cbind(1:model$n_states,max.col(i,ties.method="first")))
      paramEM<-lapply(1:model$n_channels,function(i) {
        x<-rbind(emissNZ[[i]],maxEM[[i]])
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),,drop = FALSE]
      })
      npEM<-sapply(paramEM,nrow)
    } else {
      maxEM <- lapply(model$emission_matrix,function(i) max.col(i,ties.method="first"))
      paramEM<-lapply(1:model$n_channels,function(i) {
        x<-rbind(emissNZ[[i]],c(1,maxEM[[i]]))
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE))][2]
      })
      npEM<-length(unlist(paramEM))
    }
    
    maxEMvalue<-lapply(1:model$n_channels, function(i) 
      apply(model$emission_matrix[[i]],1,max))

    
    emissNZ<-array(0,c(model$n_states,max(model$n_symbols),model$n_channels))
    for(i in 1:model$n_channels){
      emissNZ[,1:model$n_symbols[i],i]<-model$emission_matrix[[i]] > 0
      emissNZ[,1:model$n_symbols[i],i][maxEM[[i]]]<-0
      
    }       
    
    initialvalues<-c(if((npTM+sum(npEM)+npIPAll)>0) log(c(
      if(npTM>0) model$transition_matrix[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$n_channels,
        function(x) model$emission_matrix[[x]][paramEM[[x]]])),
      if(npIPAll>0) unlist(sapply(1:original_model$n_clusters,function(m)
        if(npIP[m]>0) original_model$initial_probs[[m]][paramIP[[m]]]))
    )),
      model$coefficients[,-1]
    )         
    
    
    objectivef<-function(pars,model, estimate = TRUE){      
      
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
            emissionArray[,1:model$n_symbols[i],i]/rowSums(emissionArray[,1:model$n_symbols[i],i])
        }
      }
      for(m in 1:original_model$n_clusters){
        if(npIP[m]>0){
          original_model$initial_probs[[m]][maxIP[[m]]] <- maxIPvalue[[m]] # Not needed?
          original_model$initial_probs[[m]][paramIP[[m]]] <- exp(pars[npTM+sum(npEM)+c(0,cumsum(npIP))[m]+
              1:npIP[m]])
          original_model$initial_probs[[m]][] <- original_model$initial_probs[[m]]/sum(original_model$initial_probs[[m]])
        }
      }
      model$initial_probs <- unlist(original_model$initial_probs)
      model$coefficients[,-1] <- pars[npTM+sum(npEM)+npIPAll+1:npCoef]
     
      if(estimate){
        objectivex(model$transition_matrix, emissionArray, model$initial_probs, obsArray, 
          transNZ, emissNZ, initNZ, model$n_symbols, 
          model$coefficients, model$X, model$n_states_in_clusters)
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
        lb <- c(rep(-10,length(initialvalues)-npCoef),rep(-150/apply(abs(model$X),2,max),model$n_clusters-1))
      }
      lb <- pmin(lb, 2*initialvalues)
      if(missing(ub)){
        ub <- c(rep(10,length(initialvalues)-npCoef),rep(150/apply(abs(model$X),2,max),model$n_clusters-1))
        #pmin(c(rep(250,length(initialvalues)-npCoef),rep(250/apply(abs(model$X),2,max),model$n_clusters-1)),
        #  pmax(250, 2*initialvalues))
      }
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
      ub <- c(rep(300,length(initialvalues)-npCoef),rep(300/apply(abs(model$X),2,max),model$n_clusters-1))
      ub <- pmax(ub, 2*initialvalues)
     localres<-nloptr(x0 = initialvalues, eval_f = objectivef,
        opts = control_local, model = model, estimate = TRUE, ub = ub, ...)
      
      model <- objectivef(localres$solution,model, FALSE)
      ll <- -localres$objective
    } else localres <- NULL
    
    rownames(model$coefficients) <- colnames(model$X)
    colnames(model$coefficients) <- model$cluster_names
    
  } else globalres <- localres <- NULL
  
  if(model$n_channels == 1){
    model$observations <- model$observations[[1]]
    model$emission_matrix <- model$emission_matrix[[1]]
  }

  model <- spread_models(model)
  attr(model, "df") <- df
  attr(model, "nobs") <- nobs
  
  for(i in 1:model$n_clusters){
    dimnames(model$transition_matrix[[i]]) <- dimnames(original_model$transition_matrix[[i]])
    for(j in 1:model$n_channels){
      dimnames(model$emission_matrix[[i]][[j]]) <- dimnames(original_model$emission_matrix[[i]][[j]])
    }
  }
  
  list(model = model, 
    logLik = ll, em_results=resEM[5:7], global_results = globalres, local_results = localres)
  
}