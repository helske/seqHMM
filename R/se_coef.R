#' Standard Errors for Regression Coefficients of Mixture Hidden Markov Model
#'
#' @importFrom numDeriv jacobian
#' @param model Object of class \code{mhmm}.
#' @param conditional If \code{TRUE} (default), compute standard errors using 
#' analytical formulas by assuming the other model parameters as fixed.
#' If \code{FALSE}, first the numerical approximation of the full Hessian is computed using analytical gradients, 
#' which is then inverted, and the standard errors of corresponding coefficients are extracted. 
#' Note that computing the non-conditional standard errors can be slow for large models.
#' @param ... Additional arguments to function \code{jacobian} of \code{numDeriv} package.
#' @return Matrix containing the standard errors for coefficients.
#' @export
#'
se_coef <- function(model, conditional = TRUE, ...){
  if (conditional) {
    matrix(c(rep(0,model$n_covariates),
      sqrt(diag(varcoef(model$coefficients, model$X, model$n_states)))),
      nrow = model$n_covariates, ncol = model$n_clusters)
  } else {
    # copied from fit_mhmm
    # 
    cmodel <- combine_models(model)
    
    if(cmodel$n_channels == 1){
      cmodel$observations <- list(cmodel$observations)
      cmodel$emission_matrix <- list(cmodel$emission_matrix)
    }
    
    obsArray<-array(0, c(cmodel$n_sequences, cmodel$length_of_sequences, 
      cmodel$n_channels))
    for(i in 1:cmodel$n_channels){
      obsArray[,,i]<-data.matrix(cmodel$observations[[i]])-1
      obsArray[,,i][obsArray[,,i]>cmodel$n_symbols[i]] <- cmodel$n_symbols[i]
    } 
    emissionArray<-array(1,c(cmodel$n_states,max(cmodel$n_symbols)+1,cmodel$n_channels))
    for(i in 1:cmodel$n_channels)
      emissionArray[,1:cmodel$n_symbols[i],i]<-cmodel$emission_matrix[[i]]
    
    maxIP <- maxIPvalue <- npIP <- numeric(model$n_clusters)  
    paramIP <-  initNZ <-vector("list",model$n_clusters)
    for(m in 1:model$n_clusters){
      # Index of largest initial probability
      maxIP[m] <- which.max(model$initial_probs[[m]])
      # Value of largest initial probability
      maxIPvalue[m] <- model$initial_probs[[m]][maxIP[m]]
      # Rest of non-zero probs
      paramIP[[m]] <- setdiff(which(model$initial_probs[[m]]>0),maxIP[m])
      npIP[m] <- length(paramIP[[m]])
      initNZ[[m]]<-model$initial_probs[[m]]>0
      initNZ[[m]][maxIP[m]]<-0
    }
    initNZ<-unlist(initNZ)
    npIPAll <- sum(unlist(npIP))
    # Largest transition probabilities (for each row)
    x<-which(cmodel$transition_matrix>0,arr.ind=TRUE)  
    transNZ<-x[order(x[,1]),]
    maxTM<-cbind(1:cmodel$n_states,max.col(cmodel$transition_matrix,ties.method="first"))
    maxTMvalue<-apply(cmodel$transition_matrix,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM)|duplicated(paramTM,fromLast=TRUE)),,drop=FALSE]
    npTM<-nrow(paramTM)
    transNZ<-cmodel$transition_matrix>0
    transNZ[maxTM]<-0    
    
    npCoef<-length(cmodel$coefficients[,-1])
    cmodel$coefficients[,1] <- 0
    
    
    emissNZ<-lapply(cmodel$emission_matrix,function(i){
      x<-which(i>0,arr.ind=TRUE) 
      x[order(x[,1]),]
    })
    
    if(cmodel$n_states > 1){
      maxEM <- lapply(cmodel$emission_matrix,function(i) cbind(1:cmodel$n_states,max.col(i,ties.method="first")))
      paramEM<-lapply(1:cmodel$n_channels,function(i) {
        x<-rbind(emissNZ[[i]],maxEM[[i]])
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE)),,drop = FALSE]
      })
      npEM<-sapply(paramEM,nrow)
    } else {
      maxEM <- lapply(cmodel$emission_matrix,function(i) max.col(i,ties.method="first"))
      paramEM<-lapply(1:cmodel$n_channels,function(i) {
        x<-rbind(emissNZ[[i]],c(1,maxEM[[i]]))
        x[!(duplicated(x)|duplicated(x,fromLast=TRUE))][2]
      })
      npEM<-length(unlist(paramEM))
    }
    
    maxEMvalue<-lapply(1:cmodel$n_channels, function(i) 
      apply(cmodel$emission_matrix[[i]],1,max))
    
    
    emissNZ<-array(0,c(cmodel$n_states,max(cmodel$n_symbols),cmodel$n_channels))
    for(i in 1:cmodel$n_channels){
      emissNZ[,1:cmodel$n_symbols[i],i]<-cmodel$emission_matrix[[i]] > 0
      emissNZ[,1:cmodel$n_symbols[i],i][maxEM[[i]]]<-0
      
    }       
    
    initialvalues<-c(if((npTM+sum(npEM)+npIPAll)>0) log(c(
      if(npTM>0) cmodel$transition_matrix[paramTM],
      if(sum(npEM)>0) unlist(sapply(1:cmodel$n_channels,
        function(x) cmodel$emission_matrix[[x]][paramEM[[x]]])),
      if(npIPAll>0) unlist(sapply(1:model$n_clusters,function(m)
        if(npIP[m]>0) model$initial_probs[[m]][paramIP[[m]]]))
    )),
      cmodel$coefficients[,-1]
    )         
    
    coef_ind <- npTM+sum(npEM)+npIPAll+1:npCoef
    objectivef<-function(pars,cmodel){      
      
      if(npTM>0){
        cmodel$transition_matrix[maxTM]<-maxTMvalue     
        cmodel$transition_matrix[paramTM]<-exp(pars[1:npTM])
        cmodel$transition_matrix<-cmodel$transition_matrix/rowSums(cmodel$transition_matrix)    
      }
      if(sum(npEM)>0){            
        for(i in 1:cmodel$n_channels){
          emissionArray[,1:cmodel$n_symbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]    
          emissionArray[,1:cmodel$n_symbols[i],i][paramEM[[i]]]<-
            exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])
          emissionArray[,1:cmodel$n_symbols[i],i]<-
            emissionArray[,1:cmodel$n_symbols[i],i]/rowSums(emissionArray[,1:cmodel$n_symbols[i],i])
        }
      }
      for(m in 1:model$n_clusters){
        if(npIP[m]>0){
          model$initial_probs[[m]][maxIP[[m]]] <- maxIPvalue[[m]] # Not needed?
          model$initial_probs[[m]][paramIP[[m]]] <- exp(pars[npTM+sum(npEM)+c(0,cumsum(npIP))[m]+
              1:npIP[m]])
          model$initial_probs[[m]][] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
        }
      }
      cmodel$initial_probs <- unlist(model$initial_probs)
      cmodel$coefficients[,-1] <- pars[coef_ind]
      
      
      objectivex(cmodel$transition_matrix, emissionArray, cmodel$initial_probs, obsArray, 
        transNZ, emissNZ, initNZ, cmodel$n_symbols, 
        cmodel$coefficients, cmodel$X, cmodel$n_states_in_clusters)$gradient
      
    }
    sqrt(diag(solve(jacobian(objectivef,initialvalues, ...))[coef_ind,coef_ind]))
  }
}