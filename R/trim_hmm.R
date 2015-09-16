#' Trim Small Probabilities of Hidden Markov Model
#' 
#' Function \code{trim_hmm} tries to set small insignificant probabilities to zero 
#' without decreasing the likelihood significantly.
#' 
#' @export
#' @param model Model of class \code{hmm} or \code{mhmm} for which 
#'   trimming is performed.
#' @param maxit Number of iterations. After zeroing small values, the model is 
#'   refitted, and this is repeated until there is nothing to trim or maxit 
#'   iterations are done.
#' @param return_loglik Return the log-likelihood of the trimmed model together with
#'   the model object. The default is \code{FALSE}.
#' @param zerotol Values smaller than this are trimmed to zero.
#' @param ... Further parameters passed on to \code{fit_hmm}.
#'   
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and fitting 
#' hidden Markov models; and \code{\link{hmm_biofam}} for information on the model used 
#' in the example.
#'   
#' @examples 
#' data(hmm_biofam)
#' 
#' # Testing if changing parameter values smaller than 1e-04 to zero 
#' # leads to improved log-likelihood.
#' hmm_trim <- trim_hmm(hmm_biofam, zerotol=1e-04, maxit=10)
#' 
trim_hmm <- function(model, maxit = 0, return_loglik=FALSE, zerotol=1e-8, ...){
  
  ll_original <- logLik(model)
  model_original <- model
  
  if(inherits(model, "hmm")){
    
    if(model$n_channels==1){
      if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) || 
             any(model$transition_matrix < zerotol & model$transition_matrix > 0)
           || any(model$emission_matrix < zerotol & model$emission_matrix > 0))){
        
        print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }
      
      model$initial_probs[model$initial_probs < zerotol] <- 0
      model$initial_probs <- model$initial_probs/sum(model$initial_probs)
      model$transition_matrix[model$transition_matrix < zerotol] <- 0
      model$transition_matrix <- model$transition_matrix/rowSums(model$transition_matrix)
      model$emission_matrix[model$emission_matrix < zerotol] <- 0
      model$emission_matrix <- model$emission_matrix/rowSums(model$emission_matrix)
      
      
      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_hmm(model,...)
          ll <- fit$logLik
          
          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break
          
          if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) || 
                 any(model$transition_matrix < zerotol & model$transition_matrix > 0)
               || any(model$emission_matrix < zerotol & model$emission_matrix > 0)))
            break
          
          model$initial_probs[model$initial_probs < zerotol] <- 0
          model$initial_probs <- model$initial_probs/sum(model$initial_probs)
          model$transition_matrix[model$transition_matrix < zerotol] <- 0
          model$transition_matrix <- model$transition_matrix/rowSums(model$transition_matrix)
          model$emission_matrix[model$emission_matrix < zerotol] <- 0
          model$emission_matrix <- model$emission_matrix/rowSums(model$emission_matrix)
          
        }
      }
      
    } else {
      if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) || 
             any(model$transition_matrix < zerotol & model$transition_matrix > 0)
           || any(sapply(model$emission_matrix,function(x) any(x < zerotol & x > 0))))){
        print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }
      model$initial_probs[model$initial_probs < zerotol] <- 0
      model$initial_probs <- model$initial_probs/sum(model$initial_probs)
      model$transition_matrix[model$transition_matrix < zerotol] <- 0
      model$transition_matrix <- model$transition_matrix/rowSums(model$transition_matrix)
      for(i in 1:model$n_channels){
        model$emission_matrix[[i]][model$emission_matrix[[i]] < zerotol] <- 0
        model$emission_matrix[[i]] <- model$emission_matrix[[i]]/
          rowSums(model$emission_matrix[[i]])
        
      }
      
      
      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_hmm(model,...)
          ll <- fit$logLik
          
          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break
          
          if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) || 
                 any(model$transition_matrix < zerotol & model$transition_matrix > 0)
               || any(sapply(model$emission_matrix,function(x) any(x < zerotol & x > 0)))))
            break
          
          model$initial_probs[model$initial_probs < zerotol] <- 0
          model$initial_probs <- model$initial_probs/sum(model$initial_probs)
          model$transition_matrix[model$transition_matrix < zerotol] <- 0
          model$transition_matrix <- model$transition_matrix/rowSums(model$transition_matrix)
          for(i in 1:model$n_channels){
            model$emission_matrix[[i]][model$emission_matrix[[i]] < zerotol] <- 0
            model$emission_matrix[[i]] <- model$emission_matrix[[i]]/
              rowSums(model$emission_matrix[[i]])          
          }       
        }
      }
    }
    
  }else if(inherits(model, "mhmm")){
    if(model$n_channels==1){
      if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) || 
             any(unlist(model$transition_matrix) < zerotol & unlist(model$transition_matrix) > 0)
           || any(unlist(model$emission_matrix) < zerotol & unlist(model$emission_matrix) > 0))){
        
        print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }

      for(m in 1:model$n_clusters){
        model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
        model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
        model$transition_matrix[[m]][model$transition_matrix[[m]] < zerotol] <- 0
        model$transition_matrix[[m]] <- model$transition_matrix[[m]]/rowSums(model$transition_matrix[[m]])
        model$emission_matrix[[m]][model$emission_matrix[[m]] < zerotol] <- 0
        model$emission_matrix[[m]] <- model$emission_matrix[[m]]/rowSums(model$emission_matrix[[m]])
      }
      
      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_mhmm(model,...)
          ll <- fit$logLik
          
          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break
          
          if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) || 
                 any(unlist(model$transition_matrix) < zerotol & unlist(model$transition_matrix) > 0)
               || any(unlist(model$emission_matrix) < zerotol & unlist(model$emission_matrix) > 0)))
            break
          
          for(m in 1:model$n_clusters){
            model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
            model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
            model$transition_matrix[[m]][model$transition_matrix[[m]] < zerotol] <- 0
            model$transition_matrix[[m]] <- model$transition_matrix[[m]]/rowSums(model$transition_matrix[[m]])
            model$emission_matrix[[m]][model$emission_matrix[[m]] < zerotol] <- 0
            model$emission_matrix[[m]] <- model$emission_matrix[[m]]/rowSums(model$emission_matrix[[m]])
          }  
        }
      }
      
    } else {
      if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) || 
             any(unlist(model$transition_matrix) < zerotol & unlist(model$transition_matrix) > 0)
           || any(unlist(model$emission_matrix) < zerotol & unlist(model$emission_matrix) > 0))){
        print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }
      for(m in 1:model$n_clusters){
        model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
        model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
        model$transition_matrix[[m]][model$transition_matrix[[m]] < zerotol] <- 0
        model$transition_matrix[[m]] <- model$transition_matrix[[m]]/rowSums(model$transition_matrix[[m]])
        for(i in 1:model$n_channels){
          model$emission_matrix[[m]][[i]][model$emission_matrix[[m]][[i]] < zerotol] <- 0
          model$emission_matrix[[m]][[i]] <- model$emission_matrix[[m]][[i]]/
            rowSums(model$emission_matrix[[m]][[i]])
          
        }
      }
      
      
      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_mhmm(model,...)
          ll <- fit$logLik
          
          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break
          
          if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) || 
                 any(unlist(model$transition_matrix) < zerotol & unlist(model$transition_matrix) > 0)
               || any(unlist(model$emission_matrix) < zerotol & unlist(model$emission_matrix) > 0)))
            break
          
          for(m in 1:model$n_clusters){
            model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
            model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
            model$transition_matrix[[m]][model$transition_matrix[[m]] < zerotol] <- 0
            model$transition_matrix[[m]] <- model$transition_matrix[[m]]/rowSums(model$transition_matrix[[m]])
            for(i in 1:model$n_channels){
              model$emission_matrix[[m]][[i]][model$emission_matrix[[m]][[i]] < zerotol] <- 0
              model$emission_matrix[[m]][[i]] <- model$emission_matrix[[m]][[i]]/
                rowSums(model$emission_matrix[[m]][[i]])
              
            }
          }
          
        }
      }
    }
  }else{
    stop("An object of class hmm or mhmm required.")
  }
  
  
  if(maxit > 0)
    print(paste(ii,"iteration(s) used."))
  
  if(ll0 < ll_original){
    print(paste("Log-likelihood of the trimmed model is smaller than the original log-likelihood, ll_trim-ll_orig =", signif(ll0-ll_original, 3)))
  } else print(paste("Trimming improved log-likelihood, ll_trim-ll_orig =", signif(ll0-ll_original, 3)))
  
  
  
  if(return_loglik){
    list(model=model,loglik=ll0)
  } else model 
}