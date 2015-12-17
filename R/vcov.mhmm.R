#' Variance-Covariance Matrix for Coefficients of Covariates of Mixture Hidden Markov Model
#'
#' Returns the asymptotic covariances matrix of maximum likelihood estimates of
#' the coefficients corresponding to the explanatory variables of the model.
#'
#' @details The conditional standard errors are computed using
#' analytical formulas by assuming that the coefficient estimates are not correlated with
#' other model parameter estimates (or that the other parameters are assumed to be fixed).
#' This often underestimates the true standard errors, but is substantially
#' faster approach for preliminary analysis. The non-conditional standard errors
#' are based on the numerical approximation of the full Hessian of the coefficients
#' and the model parameters corresponding to nonzero probabilities.
#' Computing the non-conditional standard errors can be slow for large models as
#' the Jacobian of analytical gradients is computed using finite difference approximation.
#'
#' @importFrom numDeriv jacobian
#' @param object Object of class \code{mhmm}.
#' @param conditional If \code{TRUE} (default), the standard errors are
#' computed conditional on other model parameters. See details.
#' @param threads Number of threads to use in parallel computing. Default is 1.
#' @param log_space Make computations using log-space instead of scaling for greater
#' numerical stability at cost of decreased computational performance. Default is \code{FALSE}.
#' @param ... Additional arguments to function \code{jacobian} of \code{numDeriv} package.
#' @return Matrix containing the variance-covariance matrix of coefficients.
#' @export
#'
vcov.mhmm <- function(object, conditional = TRUE, threads = 1, log_space = FALSE, ...){

  if (conditional) {
    vcovm <- varcoef(object$coefficients, object$X)
  } else {
    if (threads < 1) stop ("Argument threads must be a positive integer.")
    # copied from fit_model
    #
    original_model <- object
    model <- combine_models(object)

    if(model$n_channels ==  1){
      model$observations <- list(model$observations)
      model$emission_probs <- list(model$emission_probs)
    }

    obsArray <- array(0, c(model$n_sequences, model$length_of_sequences,
      model$n_channels))
    for(i in 1:model$n_channels){
      obsArray[,,i] <- data.matrix(model$observations[[i]])-1
      obsArray[,,i][obsArray[,,i]>model$n_symbols[i]] <- model$n_symbols[i]
    }
    obsArray <- aperm(obsArray)

    emissionArray <- array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
    for(i in 1:model$n_channels)
      emissionArray[,1:model$n_symbols[i],i] <- model$emission_probs[[i]]

    maxIP <- maxIPvalue <- npIP <- numeric(original_model$n_clusters)
    paramIP <- initNZ <- vector("list",original_model$n_clusters)
    for(m in 1:original_model$n_clusters){
      # Index of largest initial probability
      maxIP[m] <- which.max(original_model$initial_probs[[m]])
      # Value of largest initial probability
      maxIPvalue[m] <- original_model$initial_probs[[m]][maxIP[m]]
      # Rest of non-zero probs
      paramIP[[m]] <- setdiff(which(original_model$initial_probs[[m]]>0),maxIP[m])
      npIP[m] <- length(paramIP[[m]])
      initNZ[[m]] <- original_model$initial_probs[[m]]>0
      initNZ[[m]][maxIP[m]] <- 0
    }
    initNZ <- unlist(initNZ)
    npIPAll <- sum(unlist(npIP))
    # Largest transition probabilities (for each row)
    x <- which(model$transition_probs>0,arr.ind = TRUE)
    transNZ <- x[order(x[,1]),]
    maxTM <- cbind(1:model$n_states,max.col(model$transition_probs,ties.method = "first"))
    maxTMvalue <- apply(model$transition_probs,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast = TRUE)),,drop = FALSE]
    npTM <- nrow(paramTM)
    transNZ <- model$transition_probs>0
    transNZ[maxTM] <- 0

    npCoef <- length(model$coefficients[,-1])
    model$coefficients[,1] <- 0


    emissNZ <- lapply(model$emission_probs,function(i){
      x <- which(i>0,arr.ind = TRUE)
      x[order(x[,1]),]
    })

    if(model$n_states > 1){
      maxEM <- lapply(model$emission_probs,function(i) cbind(1:model$n_states,max.col(i,ties.method = "first")))
      paramEM <- lapply(1:model$n_channels,function(i) {
        x <- rbind(emissNZ[[i]],maxEM[[i]])
        x[!(duplicated(x)|duplicated(x,fromLast = TRUE)),,drop = FALSE]
      })
      npEM <- sapply(paramEM,nrow)
    } else {
      maxEM <- lapply(model$emission_probs,function(i) max.col(i,ties.method = "first"))
      paramEM <- lapply(1:model$n_channels,function(i) {
        x <- rbind(emissNZ[[i]],c(1,maxEM[[i]]))
        x[!(duplicated(x)|duplicated(x,fromLast = TRUE))][2]
      })
      npEM <- length(unlist(paramEM))
    }

    maxEMvalue <- lapply(1:model$n_channels, function(i)
      apply(model$emission_probs[[i]],1,max))


    emissNZ <- array(0,c(model$n_states,max(model$n_symbols),model$n_channels))
    for(i in 1:model$n_channels){
      emissNZ[,1:model$n_symbols[i],i] <- model$emission_probs[[i]] > 0
      emissNZ[,1:model$n_symbols[i],i][maxEM[[i]]] <- 0

    }

    initialvalues <- c(if((npTM+sum(npEM)+npIPAll)>0) log(c(
      if(npTM>0) model$transition_probs[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:model$n_channels,
        function(x) model$emission_probs[[x]][paramEM[[x]]])),
      if(npIPAll>0) unlist(sapply(1:original_model$n_clusters,function(m)
        if(npIP[m]>0) original_model$initial_probs[[m]][paramIP[[m]]]))
    )),
      model$coefficients[,-1]
    )

    coef_ind <- npTM+sum(npEM)+npIPAll+1:npCoef

    objectivef <- function(pars,model){

      if(npTM>0){
        model$transition_probs[maxTM] <- maxTMvalue
        model$transition_probs[paramTM] <- exp(pars[1:npTM])
        model$transition_probs <- model$transition_probs/rowSums(model$transition_probs)
      }
      if(sum(npEM)>0){
        for(i in 1:model$n_channels){
          emissionArray[,1:model$n_symbols[i],i][maxEM[[i]]] <- maxEMvalue[[i]]
          emissionArray[,1:model$n_symbols[i],i][paramEM[[i]]] <- 
            exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])
          emissionArray[,1:model$n_symbols[i],i] <- 
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
      model$coefficients[,-1] <- pars[coef_ind]

      if (!log_space) {
        objectivex(model$transition_probs, emissionArray, model$initial_probs, obsArray,
          transNZ, emissNZ, initNZ, model$n_symbols,
          model$coefficients, model$X, model$n_states_in_clusters, threads)$gradient
      } else {
        log_objectivex(model$transition_probs, emissionArray, model$initial_probs, obsArray,
          transNZ, emissNZ, initNZ, model$n_symbols,
          model$coefficients, model$X, model$n_states_in_clusters, threads)$gradient
      }
    }
    vcovm <- solve(jacobian(objectivef, initialvalues, model = model, ...))[coef_ind, coef_ind]
  }
  rownames(vcovm) <- colnames(vcovm) <- paste(
    rep(object$cluster_names[-1], each = object$n_covariates),
    rownames(object$coefficients), sep = ": ")
  vcovm
}
