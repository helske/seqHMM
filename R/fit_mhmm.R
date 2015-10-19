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
#'   \code{\link{nloptr}} (possibly after the EM step). The default is \code{FALSE}.
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
#'   }
#' @param threads Number of threads to use in parallel computing. Default is 1.
#' @param ... Additional arguments to nloptr
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
#'   a local optimization (LBFGS by default) is run with a lower tolerance to find the 
#'   optimum with high precision.
#'   
#'   There are some theoretical guarantees that the MLSL method used as the default 
#'   optimizer in step 2 shoud find all local optima in a finite number of local 
#'   optimizations. Of course, it might not always succeed in a reasonable time. 
#'   The EM algorithm can help in finding good boundaries for the search, especially 
#'   with good starting values, but in some cases it can mislead. A good strategy is to 
#'   try a couple of different fitting options with different combinations of the methods: 
#'   e.g. all steps, only global and local steps, and a few evaluations of EM followed by 
#'   global and local optimization.
#'   
#'   By default, the estimation time is limited to 60 seconds in steps 2 and 3, so it is 
#'   advisable to change the default settings for the final analysis. 
#'   
#'   Any method available in the \code{nloptr} function can be used for the global and 
#'   local steps.
#' @return List with components \item{model}{Estimated model. } 
#'   \item{logLik}{Log-likelihood of the estimated model. } 
#'   \item{em_results}{Results after the EM step. } 
#'   \item{global_results}{Results after the global step. }
#'   \item{local_results}{Results after the local step. }
#' @seealso \code{\link{build_mhmm}} for building MHMMs; \code{\link{summary.mhmm}}
#'   for a summary of a MHMM; \code{\link{separate_mhmm}} for reorganizing a MHMM into 
#'   a list of separate hidden Markov models;
#'   \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and
#'   fitting hidden Markov models; \code{\link{plot.mhmm}}
#'   for plotting \code{mhmm} objects; and \code{\link{mssplot}} for plotting
#'   stacked sequence plots of \code{mhmm} objects.
#' @examples
#' 
#' # Single-channel
#' 
#' data(mvad, package = "TraMineR")
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
#' init_mhmm_mvad <- build_mhmm(
#'   observations = mvad.seq, 
#'   transition_probs = list(tr_1, tr_2), 
#'   emission_probs = list(emiss_1, emiss_2), 
#'   initial_probs = list(init_1, init_2))                   
#' 
#' # Fitting the model with the EM algorithm
#' fit_mvad <- fit_mhmm(
#'   init_mhmm_mvad, local_step = FALSE)
#' fit_mvad$logLik # -20639.1
#' \dontrun{
#' # Run EM algorithm 25 times with simulated starting values
#' set.seed(321)
#' fit_mvad2 <- fit_mhmm(init_mhmm_mvad, control_em=list(restarts = 25))
#' fit_mvad2$logLik # -14651.43
#' }
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
#' init_mhmm_bf <- build_mhmm(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_probs = list(trans_1, trans_1, trans_2),
#'   emission_probs = list(list(emiss_1_child, emiss_1_marr, emiss_1_left),
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
#'   init_mhmm_bf, em_step = TRUE, global_step = FALSE, local_step = FALSE)
#' mhmm_1$logLik # -3081.383
#' 
#' \dontrun{
#' # EM with LBFGS
#' mhmm_2 <- fit_mhmm(
#'   init_mhmm_bf, em_step = TRUE, global_step = FALSE, local_step = TRUE)
#' mhmm_2$logLik # -3081.383
#' 
#' # Only LBFGS
#' mhmm_3 <- fit_mhmm(
#'   init_mhmm_bf, em_step = FALSE, global_step = FALSE, local_step = TRUE,
#'   control_local = list(maxeval = 5000, maxtime = 0))
#' mhmm_3$logLik # -3087.499373
#' 
#' # Global optimization via MLSL_LDS with LBFGS as local optimizer and final polisher
#' mhmm_4 <- fit_mhmm(
#'   init_mhmm_bf, em_step = FALSE, global_step = TRUE, local_step = TRUE, 
#'   control_global = list(maxeval = 5000, maxtime = 0))
#' mhmm_4$logLik # -3150.796
#' 
#' # As previously, but now we use ten iterations from EM algorithm for 
#' # defining initial values and boundaries
#' # Note smaller maxeval for global optimization
#' mhmm_5 <- fit_mhmm(
#'   init_mhmm_bf, em_step = TRUE, global_step = TRUE, local_step = TRUE, 
#'   control_em = list(maxeval = 10), control_global = list(maxeval = 1000, maxtime = 0),
#'   control_local = list(maxeval = 500, maxtime = 0))
#' mhmm_5$logLik #-3081.383
#' }
#' 

fit_mhmm <- function(model, em_step = TRUE, global_step = FALSE, local_step = TRUE, 
  control_em = list(), control_global = list(), control_local = list(), lb, ub, threads = 1, ...){

  
  if(!inherits(model, "mhmm"))
    stop("Argument model must be an object of class 'mhmm'.")
  
  if(!em_step && !global_step && !local_step){
    stop("No method chosen for estimation. Choose at least one from em_step, global_step, and local_step.")
  }
  if (threads < 1) stop ("Argument threads must be a positive integer.")
  
  df <- attr(model, "df")
  nobs <- attr(model, "nobs")
  original_model <- model
  model <- combine_models(model)
  
  if(model$n_channels == 1){
    model$observations <- list(model$observations)
    model$emission_probs <- list(model$emission_probs)
  }
  
  obsArray<-array(0, c(model$n_sequences, model$length_of_sequences, 
    model$n_channels))
  for(i in 1:model$n_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$n_symbols[i]] <- model$n_symbols[i]
  } 
  emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
  for(i in 1:model$n_channels)
    emissionArray[,1:model$n_symbols[i],i]<-model$emission_probs[[i]]
  
  if(em_step){
    em.con <- list(print_level = 0, maxeval=100, reltol=1e-8, restarts = 0,
      restart_transition = TRUE, restart_emission = TRUE, sd_restart = 0.25)
    nmsC <- names(em.con)  
    em.con[(namc <- names(control_em))] <- control_em
    if (length(noNms <- namc[!namc %in% nmsC])) 
      warning("Unknown names in control_em: ", paste(noNms, collapse = ", "))
    
    resEM <- EMx(model$transition_probs, emissionArray, model$initial_probs, obsArray, 
      model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters, 
      em.con$maxeval, em.con$reltol,em.con$print_level, threads)
    
    if(resEM$error != 0){
      err_msg <- switch(resEM$error, 
        "1" = "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
        "2" = "Estimation of coefficients of covariates failed due to singular Hessian.",
        "3" = "Estimation of coefficients of covariates failed due to non-finite cluster probabilities.",
        "4" = "Non-finite log-likelihood.")
      if(!global_step && !local_step && em.con$restarts == 0){
        stop(paste("EM algorithm failed:", err_msg))
      } else warning(paste("EM algorithm failed:", err_msg))
      resEM$logLik <- -Inf
    }
    
    if (em.con$restarts > 0 & (em.con$restart_transition | em.con$restart_emission)) {
      if (resEM$error == 0) {
        random_emiss <- resEM$emissionArray
        random_emiss[(random_emiss < 1e-4) & (emissionArray >= 1e-4)] <- 1e-4
        for (j in 1:model$n_channels) {
          random_emiss[,1:model$n_symbols[j],j] <- random_emiss[,1:model$n_symbols[j],j] / rowSums(random_emiss[,1:model$n_symbols[j],j])
        }
        random_trans <- resEM$transitionMatrix
        random_trans[(random_trans < 1e-4) & (model$transition_probs >= 1e-4)] <- 1e-4
        random_trans <- random_trans / rowSums(random_trans)
      } else {
        random_emiss <- emissionArray
        random_trans <- model$transition_probs
      }
      
      if (em.con$restart_transition) {
        nz_trans <- (random_trans > 0 & random_trans < 1)
        np_trans <- sum(nz_trans)
        base_trans <- random_trans[nz_trans]
      }
      if (em.con$restart_emission) {
        nz_emiss <- (random_emiss > 0 & random_emiss < 1)
        np_emiss <- sum(nz_emiss)
        base_emiss <- random_emiss[nz_emiss]
      }
      
      for (i in 1:em.con$restarts) {
        if (em.con$restart_transition) {
          random_trans[nz_trans] <- abs(base_trans + rnorm(np_trans, sd = em.con$sd_restart))
          random_trans <- random_trans / rowSums(random_trans)
        }
        if (em.con$restart_emission) {
          random_emiss[nz_emiss] <- abs(base_emiss + rnorm(np_emiss, sd = em.con$sd_restart))
          for (j in 1:model$n_channels) {
            random_emiss[,1:model$n_symbols[j],j] <- random_emiss[,1:model$n_symbols[j],j] / rowSums(random_emiss[,1:model$n_symbols[j],j])
          }
        }
        resEMi <- EMx(random_trans, random_emiss, model$initial_probs, obsArray, 
          model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters, em.con$maxeval, 
          em.con$reltol,em.con$print_level, threads)
        if (resEMi$error != 0) {
          err_msg <- switch(resEMi$error, 
            "1" = "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
            "2" = "Estimation of coefficients of covariates failed due to singular Hessian.",
            "3" = "Estimation of coefficients of covariates failed due to non-finite cluster probabilities.",
            "4" = "Non-finite log-likelihood.")
          warning(paste("EM algorithm failed:", err_msg))
        } else {
          if (resEMi$logLik > resEM$logLik) {
            resEM <- resEMi
          }
        }
      }
      
    }
    if (resEM$error == 0) {
      if (resEM$change < -1e-5)
        warning("EM algorithm stopped due to decreasing log-likelihood. ")
      
      emissionArray <- resEM$emissionArray
      for (i in 1:model$n_channels)
        model$emission_probs[[i]][] <- emissionArray[ , 1:model$n_symbols[i], i]
      
      if(global_step || local_step){
        k <- 0
        for(m in 1:model$n_clusters){
          original_model$initial_probs[[m]] <- unname(resEM$initialProbs[(k+1):(k+model$n_states_in_clusters[m])])
          k <- sum(model$n_states_in_clusters[1:m])
        }
      } else {
        model$initial_probs <- resEM$initialProbs
      }
      
      model$transition_probs[]<-resEM$transitionMatrix
      model$coefficients[]<-resEM$coefficients
      ll <- resEM$logLik
    } else resEM <-NULL
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
    x<-which(model$transition_probs>0,arr.ind=TRUE)  
    transNZ<-x[order(x[,1]),]
    maxTM<-cbind(1:model$n_states,max.col(model$transition_probs,ties.method="first"))
    maxTMvalue<-apply(model$transition_probs,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),,drop=FALSE]
    npTM<-nrow(paramTM)
    transNZ<-model$transition_probs>0
    transNZ[maxTM]<-0    
    
    npCoef<-length(model$coefficients[,-1])
    model$coefficients[,1] <- 0
    
    
    emissNZ<-lapply(model$emission_probs,function(i){
      x<-which(i>0,arr.ind=TRUE) 
      x[order(x[,1]),]
    })
    
    if(model$n_states > 1){
      maxEM <- lapply(model$emission_probs,function(i) cbind(1:model$n_states,max.col(i,ties.method="first")))
      paramEM<-lapply(1:model$n_channels,function(i) {
        x<-rbind(emissNZ[[i]],maxEM[[i]])
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),,drop = FALSE]
      })
      npEM<-sapply(paramEM,nrow)
    } else {
      maxEM <- lapply(model$emission_probs,function(i) max.col(i,ties.method="first"))
      paramEM<-lapply(1:model$n_channels,function(i) {
        x<-rbind(emissNZ[[i]],c(1,maxEM[[i]]))
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE))][2]
      })
      npEM<-length(unlist(paramEM))
    }
    
    maxEMvalue<-lapply(1:model$n_channels, function(i) 
      apply(model$emission_probs[[i]],1,max))
    
    
    emissNZ<-array(0,c(model$n_states,max(model$n_symbols),model$n_channels))
    for(i in 1:model$n_channels){
      emissNZ[,1:model$n_symbols[i],i]<-model$emission_probs[[i]] > 0
      emissNZ[,1:model$n_symbols[i],i][maxEM[[i]]]<-0
      
    }       
    
    initialvalues<-c(if((npTM+sum(npEM)+npIPAll)>0) log(c(
      if(npTM>0) model$transition_probs[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$n_channels,
        function(x) model$emission_probs[[x]][paramEM[[x]]])),
      if(npIPAll>0) unlist(sapply(1:original_model$n_clusters,function(m)
        if(npIP[m]>0) original_model$initial_probs[[m]][paramIP[[m]]]))
    )),
      model$coefficients[,-1]
    )         
    
    
    objectivef<-function(pars,model, estimate = TRUE){      
      
      
      if(npTM>0){
        model$transition_probs[maxTM]<-maxTMvalue     
        model$transition_probs[paramTM]<-exp(pars[1:npTM])
        model$transition_probs<-model$transition_probs/rowSums(model$transition_probs)    
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
        objectivex(model$transition_probs, emissionArray, model$initial_probs, obsArray, 
          transNZ, emissNZ, initNZ, model$n_symbols, 
          model$coefficients, model$X, model$n_states_in_clusters, threads)
      } else {
        if(sum(npEM)>0){
          for(i in 1:model$n_channels){
            model$emission_probs[[i]][]<-emissionArray[,1:model$n_symbols[i],i]
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
      }
      ub <- pmin(pmax(ub, 2*initialvalues),500)
      if(is.null(control_global$maxeval)){
        control_global$maxeval <- 10000
      }
      if(is.null(control_global$maxtime)){
        control_global$maxtime <- 60
      }
      if(is.null(control_global$algorithm)){
        control_global$algorithm <- "NLOPT_GD_MLSL_LDS"
        if(is.null(control_global$local_opts)) control_global$local_opts <- list(algorithm = "NLOPT_LD_LBFGS",  xtol_rel = 1e-4)
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
      if(is.null(control_local$algorithm)){
        control_local$algorithm <- "NLOPT_LD_LBFGS"
        if(is.null(control_local$xtol_rel)) control_local$xtol_rel <- 1e-8
      }
      ub <- c(rep(300,length(initialvalues)-npCoef),rep(300/apply(abs(model$X),2,max),model$n_clusters-1))
      ub <- pmin(pmax(ub, 2*initialvalues),500)
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
    model$emission_probs <- model$emission_probs[[1]]
  }
  
  model <- spread_models(model)
  attr(model, "df") <- df
  attr(model, "nobs") <- nobs
  
  for(i in 1:model$n_clusters){
    dimnames(model$transition_probs[[i]]) <- dimnames(original_model$transition_probs[[i]])
    for(j in 1:model$n_channels){
      dimnames(model$emission_probs[[i]][[j]]) <- dimnames(original_model$emission_probs[[i]][[j]])
    }
  }
  suppressWarnings(try(model <- trim_hmm(model, verbose = FALSE), silent = TRUE))
  list(model = model, 
    logLik = ll, em_results=resEM[5:7], global_results = globalres, local_results = localres)
  
}