#' Estimate Parameters of Hidden Markov Model
#'
#' Function \code{fit_hmm} estimates the initial state, transition and emission
#' probabilities of hidden Markov model. Initial values for estimation are taken from the
#' corresponding components of the model with preservation of original zero
#' probabilities. By default, the estimation start with EM algorithm and then switches
#' to direct numerical maximization.
#'
#' @export
#' @param model Hidden Markov model.
#' @param em_step Logical, use EM algorithm at the start of parameter estimation.
#'   The default is \code{TRUE}. Note that EM algorithm is faster than direct numerical optimization,
#'   but is even more prone to get stuck in a local optimum.
#' @param global_step Logical, use global optimization via
#'   \code{\link{nloptr}} (possibly after the EM step). The default is \code{FALSE}.
#'@param local_step Logical, use local optimization via
#'   \code{\link{nloptr}} (possibly after the EM and/or global steps). The default is \code{FALSE}.
#' @param control_em Optional list of control parameters for for EM algorithm.
#'   Possible arguments are \describe{
#'   \item{maxeval}{Maximum number of iterations, default is 1000.}
#'   \item{print_level}{Level of printing. Possible values are 0
#'   (prints nothing), 1 (prints information at start and end of algorithm), and
#'   2 (prints at every iteration).}
#'   \item{reltol}{Relative tolerance for convergence defined as \eqn{(logLik_new - logLik_old)/(abs(logLik_old)+0.1)}.
#'   Default is 1e-12.}
#'   \item{restart}{Optional list containing options for possible EM restarts with following components:
#'   \describe{
#'   \item{times}{Number of restarts of EM algorithm using random initial values. Default is 0 i.e. no restarts. }
#'   \item{transition}{Logical, should the initial transition probabilities be varied. Default is \code{TRUE}. }
#'   \item{emission}{Logical, should the initial emission probabilities be varied. Default is \code{TRUE}. }
#'   \item{sd}{Standard deviation for \code{rnorm} used in randomizing. Default is 0.25.}
#'   \item{maxeval}{Maximum number of iterations, default is 100.}
#'   \item{print_level}{Level of printing in restarted EM steps. Defaults to \code{control_em$print_level}. }
#'   \item{reltol}{Relative tolerance for convergence in restarted EM steps. Default is 1e-8.}
#'   }
#'   }
#'   }
#' @param control_global Optional list of additional arguments for
#'   \code{\link{nloptr}} argument \code{opts}. The default values are
#'   \describe{
#'    \item{algorithm}{\code{"NLOPT_GD_MLSL_LDS"}}
#'    \item{local_opts}{\code{list(algorithm = "NLOPT_LD_LBFGS",  xtol_rel = 1e-4)}}
#'    \item{maxeval}{\code{10000} (maximum number of iterations in global optimization algorithm)}
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
#'   }
#' @param threads Number of threads to use in parallel computing. Default is 1.
#' @param log_space Make computations using log-space instead of scaling for greater
#' numerical stability at cost of decreased computational performance. Default is \code{FALSE}.
#' @param ... Additional arguments to nloptr
#' @return List with components \item{model}{Estimated model. }
#'   \item{logLik}{Log-likelihood of the estimated model. }
#'   \item{em_results}{Results after the EM step. }
#'   \item{global_results}{Results after the global step. }
#'   \item{local_results}{Results after the local step. }
#' @seealso \code{\link{build_hmm}} for building Hidden Markov models before
#'   fitting, \code{\link{trim_hmm}} for finding better models by changing small
#'   parameter values to zero, and
#'   \code{\link{plot.hmm}} and \code{\link{ssplot}} for plotting
#'   hmm objects.
#' @details The fitting function provides three estimation steps: 1) EM algorithm,
#'   2) global optimization, and 3) local optimization. The user can call for one method
#'   or any combination of these steps, but should note that they are preformed in the
#'   above-mentioned order. The results from a former step are used as starting values
#'   in a latter.
#'
#'   By default the \code{fit_hmm} function starts with the EM algorithm,
#'   and finishes with LBFGS as the local optimizer.
#'
#'   It is possible to rerun EM algorithm automatically using random starting
#'   values based on the first run of EM. Number of restarts is defined by
#'   argument \code{restarts} in \code{control_em}. As EM algorithm is relatively fast, this method
#'   might be preferred option compared to proper global optimization strategy of step 2.
#'
#'   Default global optimization method (triggered via \code{global_step = TRUE} is
#'   the multilevel single-linkage method (MLSL) with the LDS modification (\code{NLOPT_GD_MLSL_LDS} as
#'   \code{algorithm} in \code{control_global}), with LBFGS as the local optimizer.
#'   The MLSL method draws random starting points and performs a local optimization
#'   from each. The LDS modification uses low-discrepancy sequences instead of
#'   pseudo-random numbers as starting points and should improve the convergence rate.
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
#'   By default, the estimation time is limited to 60 seconds in global optimization step, so it is
#'   advisable to change the default settings for the proper global optimization.
#'
#'   Any method available in the \code{nloptr} function can be used for the global and
#'   local steps.
#' @examples
#'
#' # Three-state three-channel hidden Markov model
#' # See ?hmm_biofam for five-state version
#'
#' data(biofam3c)
#'
#' # Building sequence objects
#' marr.seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child.seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left.seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' # Define colors
#' attr(marr.seq, "cpal") <- c("violetred2", "darkgoldenrod2", "darkmagenta")
#' attr(child.seq, "cpal") <- c("darkseagreen1", "coral3")
#' attr(left.seq, "cpal") <- c("lightblue", "red3")
#'
#' # Starting values for emission matrices
#'
#' emiss_marr <- matrix(NA, nrow = 3, ncol = 3)
#' emiss_marr[1,] <- seqstatf(marr.seq[, 1:5])[, 2] + 1
#' emiss_marr[2,] <- seqstatf(marr.seq[, 6:10])[, 2] + 1
#' emiss_marr[3,] <- seqstatf(marr.seq[, 11:16])[, 2] + 1
#' emiss_marr <- emiss_marr / rowSums(emiss_marr)
#'
#' emiss_child <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 1
#' emiss_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 1
#' emiss_child[3,] <- seqstatf(child.seq[, 11:16])[, 2] + 1
#' emiss_child <- emiss_child / rowSums(emiss_child)
#'
#' emiss_left <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_left[1,] <- seqstatf(left.seq[, 1:5])[, 2] + 1
#' emiss_left[2,] <- seqstatf(left.seq[, 6:10])[, 2] + 1
#' emiss_left[3,] <- seqstatf(left.seq[, 11:16])[, 2] + 1
#' emiss_left <- emiss_left / rowSums(emiss_left)
#'
#' # Starting values for transition matrix
#' trans <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#'
#' # Starting values for initial state probabilities
#' inits <- c(0.9, 0.09, 0.01)
#'
#' # Building hidden Markov model with initial parameter values
#' init_hmm_bf <- build_hmm(
#'   observations = list(marr.seq, child.seq, left.seq),
#'   transition_probs = trans,
#'   emission_probs = list(emiss_marr, emiss_child, emiss_left),
#'   initial_probs = inits)
#'
#' # Fitting the model with different optimization schemes
#'
#' # Only EM with default values
#' hmm_1 <- fit_hmm(init_hmm_bf)
#' hmm_1$logLik # -24179.1
#'
#' \dontrun{
#'
#' # Only LBFGS
#' hmm_2 <- fit_hmm(init_hmm_bf, em_step = FALSE, local_step = TRUE)
#' hmm_2$logLik # -22267.75
#'
#' # Global optimization via MLSL_LDS with LBFGS as local optimizer and final polisher
#' # This can be slow, use parallel computing by adjusting threads argument
#' # (here threads = 1 for portability issues)
#' hmm_3 <- fit_hmm(
#'   init_hmm_bf, em_step = FALSE, global_step = TRUE, local_step = TRUE,
#'   control_global = list(maxeval = 5000, maxtime = 0), threads = 1)
#' hmm_3$logLik # -22267.75
#'
#' # EM with restarts, much faster than MLSL
#' set.seed(123)
#' hmm_4 <- fit_hmm(init_hmm_bf, control_em = list(restarts = 5, print_level = 1), threads = 1)
#' hmm_4$logLik # -21675.4
#'
#' #' # Global optimization via STOGO with LBFGS as local optimizer and final polisher
#' # This can be slow, use parallel computing by adjusting threads argument
#' # (here threads = 1 for portability issues)
#' set.seed(123)
#' hmm_5 <- fit_hmm(
#'    init_hmm_bf, em_step = FALSE, global_step = TRUE, local_step = TRUE,
#' control_global = list(algorithm = "NLOPT_GD_STOGO", maxeval = 2500, maxtime = 0), threads = 1)
#' hmm_5$logLik # -22131.76
#'
#' # Using log_space results better optimum (same is also true for MLSL example above):
#' #' set.seed(123)
#' hmm_5b <- fit_hmm(
#'    init_hmm_bf, em_step = FALSE, global_step = TRUE, local_step = TRUE, log_space = TRUE,
#' control_global = list(algorithm = "NLOPT_GD_STOGO", maxeval = 2500, maxtime = 0), threads = 1)
#' hmm_5b$logLik # -21675.4

#' }
#'
fit_hmm <- function(model, em_step = TRUE, global_step = FALSE, local_step = FALSE,
  control_em=list(), control_global=list(), control_local=list(), lb, ub,
  threads = 1, log_space = FALSE, ...){

  if(!inherits(model, "hmm"))
    stop("Argument model must be an object of class 'hmm'.")

  if(!em_step && !global_step && !local_step){
    stop("No method chosen for estimation. Choose at least one from em_step, global_step, and local_step.")
  }
  if (threads < 1) stop ("Argument threads must be a positive integer.")

  if(model$n_channels == 1){
    model$observations <- list(model$observations)
    model$emission_probs <- list(model$emission_probs)
  }

  obsArray<-array(0,c(model$n_sequences,model$length_of_sequences,model$n_channels))
  for(i in 1:model$n_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$n_symbols[i]]<-model$n_symbols[i]
  }

  if (em_step) {
    em.con <- list(print_level = 0, maxeval = 1000, reltol = 1e-12, restart = NULL)
    nmsC <- names(em.con)
    em.con[(namc <- names(control_em))] <- control_em
    if (length(noNms <- namc[!namc %in% nmsC]))
      warning("Unknown names in control_em: ", paste(noNms, collapse = ", "))

    emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
    for(i in 1:model$n_channels)
      emissionArray[,1:model$n_symbols[i],i]<-model$emission_probs[[i]]
    if (!log_space) {
      resEM<-EM(model$transition_probs, emissionArray, model$initial_probs, obsArray,
        model$n_symbols, em.con$maxeval, em.con$reltol,em.con$print_level, threads)
    } else {
      resEM <- log_EM(model$transition_probs, emissionArray, model$initial_probs, obsArray,
        model$n_symbols, em.con$maxeval, em.con$reltol,em.con$print_level, threads)
    }

    if (!is.null(em.con$restart)) {
      restart.con <- list(times = 0, print_level = em.con$print_level, maxeval = 100, reltol = 1e-8,
        transition = TRUE, emission = TRUE, sd = 0.25)
      nmsC <- names(restart.con)
      restart.con[(namc <- names(control_em$restart))] <- control_em$restart
      if (length(noNms <- namc[!namc %in% nmsC]))
        warning("Unknown names in control_em$restart: ", paste(noNms, collapse = ", "))
    }


    if (!is.null(em.con$restart) && restart.con$times > 0 &&
        (restart.con$transition | restart.con$emission)) {
      random_emiss <- resEM$emissionArray
      random_emiss[(random_emiss < 1e-4) & (emissionArray >= 1e-4)] <- 1e-4
      for (j in 1:model$n_channels) {
        random_emiss[,1:model$n_symbols[j],j] <- random_emiss[,1:model$n_symbols[j],j] / rowSums(random_emiss[,1:model$n_symbols[j],j])
      }
      random_trans <- resEM$transitionMatrix
      random_trans[(random_trans < 1e-4) & (model$transition_probs >= 1e-4)] <- 1e-4
      random_trans <- random_trans / rowSums(random_trans)


      if (restart.con$transition) {
        nz_trans <- (random_trans > 0 & random_trans < 1)
        np_trans <- sum(nz_trans)
        base_trans <- random_trans[nz_trans]
      }
      if (restart.con$emission) {
        nz_emiss <- (random_emiss > 0 & random_emiss < 1)
        np_emiss <- sum(nz_emiss)
        base_emiss <- random_emiss[nz_emiss]
      }

      for (i in 1:restart.con$times) {
        if (restart.con$transition) {
          random_trans[nz_trans] <- abs(base_trans + rnorm(np_trans, sd = restart.con$sd))
          random_trans <- random_trans / rowSums(random_trans)
        }
        if (restart.con$emission) {
          random_emiss[nz_emiss] <- abs(base_emiss + rnorm(np_emiss, sd = restart.con$sd))
          for (j in 1:model$n_channels) {
            random_emiss[,1:model$n_symbols[j],j] <- random_emiss[,1:model$n_symbols[j],j] / rowSums(random_emiss[,1:model$n_symbols[j],j])
          }
        }
        if (!log_space) {
          resEMi <- EM(random_trans, random_emiss, model$initial_probs, obsArray,
            model$n_symbols, restart.con$maxeval, restart.con$reltol,restart.con$print_level, threads)
        } else {
          resEMi <- log_EM(random_trans, random_emiss, model$initial_probs, obsArray,
            model$n_symbols, restart.con$maxeval, restart.con$reltol,restart.con$print_level, threads)
        }

        if (resEMi$logLik > resEM$logLik) {
          resEM <- resEMi
        }

      }
      if (em.con$reltol < restart.con$reltol) {
        if (!log_space) {
          resEM <- EM(resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
            model$n_symbols, restart.con$maxeval, restart.con$reltol,em.con$print_level, threads)
        } else {
          resEM <- log_EM(resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
            model$n_symbols, restart.con$maxeval, restart.con$reltol,em.con$print_level, threads)
        }
      }
    }

    if (resEM$change < -1e-5)
      warning("EM algorithm stopped due to the decreasing log-likelihood. ")


    for (i in 1:model$n_channels)
      model$emission_probs[[i]][] <- resEM$emissionArray[ , 1:model$n_symbols[i], i]


    model$initial_probs[] <- resEM$initialProbs
    model$transition_probs[] <- resEM$transitionMatrix
    ll <- resEM$logLik
  } else resEM <-NULL

  if(global_step || local_step){
    maxIP<-which.max(model$initial_probs)
    maxIPvalue<-model$initial_probs[maxIP]
    paramIP<-setdiff(which(model$initial_probs>0),maxIP)

    npIP<-length(paramIP)
    initNZ<-model$initial_probs>0
    initNZ[maxIP]<-0


    x<-which(model$transition_probs>0,arr.ind=TRUE)
    transNZ<-x[order(x[,1]),]
    maxTM<-cbind(1:model$n_states,max.col(model$transition_probs,ties.method="first"))
    maxTMvalue<-apply(model$transition_probs,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),,drop=FALSE]
    npTM<-nrow(paramTM)
    transNZ<-model$transition_probs>0
    transNZ[maxTM]<-0

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
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),2]
      })
      npEM<-sapply(paramEM, length)
    }
    maxEMvalue<-lapply(1:model$n_channels, function(i)
      apply(model$emission_probs[[i]],1,max))




    emissNZ<-array(0,c(model$n_states,max(model$n_symbols),model$n_channels))
    for(i in 1:model$n_channels){
      emissNZ[,1:model$n_symbols[i],i]<-model$emission_probs[[i]] > 0
      emissNZ[,1:model$n_symbols[i],i][maxEM[[i]]]<-0
    }

    initialvalues<-c(if((npTM+sum(npEM)+npIP)>0) log(c(
      if(npTM>0) model$transition_probs[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$n_channels,
        function(x) model$emission_probs[[x]][paramEM[[x]]])),
      if(npIP>0) model$initial_probs[paramIP]))
    )

    emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
    for(i in 1:model$n_channels)
      emissionArray[,1:model$n_symbols[i],i]<-model$emission_probs[[i]]

    objectivef<-function(pars, model, estimate = TRUE){

      if(any(!is.finite(exp(pars))) && estimate)
        return(.Machine$double.xmax^075)
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
            emissionArray[,1:model$n_symbols[i],i]/rowSums(emissionArray[,1:model$n_symbols[i],i, drop = FALSE])
        }
      }

      if(npIP>0){
        model$initial_probs[maxIP]<-maxIPvalue
        model$initial_probs[paramIP]<-exp(pars[(npTM+sum(npEM)+1):(npTM+sum(npEM)+npIP)])
        model$initial_probs[]<-model$initial_probs/sum(model$initial_probs)
      }
      if(estimate){
        if (!log_space) {
          objective(model$transition_probs, emissionArray, model$initial_probs, obsArray,
            transNZ, emissNZ, initNZ, model$n_symbols, threads)
        } else {
          log_objective(model$transition_probs, emissionArray, model$initial_probs, obsArray,
            transNZ, emissNZ, initNZ, model$n_symbols, threads)
        }
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
        lb <- -50
      }
      if(missing(ub)){
        ub <- 50
      }
      lb <- pmin(lb, 2*initialvalues)
      ub <- pmin(pmax(ub, 2*initialvalues), 500)
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
        opts = control_global, model = model, estimate = TRUE,  ...)
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
        control_local$xtol_rel <- 1e-8
      }
      ub <- pmin(pmax(300, 2*initialvalues), 500)
      localres<-nloptr(x0 = initialvalues,
        eval_f = objectivef,
        opts = control_local, model = model, estimate = TRUE, ub = ub, ...)
      model <- objectivef(localres$solution,model, FALSE)
      ll <- -localres$objective
    } else localres <- NULL


  } else globalres <- localres <- NULL

  if(model$n_channels == 1){
    model$observations <- model$observations[[1]]
    model$emission_probs <- model$emission_probs[[1]]
  }

  suppressWarnings(try(model <- trim_hmm(model, verbose = FALSE), silent = TRUE))
  list(model = model, logLik = ll,
    em_results = resEM[4:6], global_results = globalres, local_results = localres)
}
