#' Trim Small Probabilities of Hidden Markov Model
#'
#' Function \code{trim_model} tries to set small insignificant probabilities to zero
#' without decreasing the likelihood.
#'
#' @export
#' @param model Model of class \code{hmm} or \code{mhmm} for which
#'   trimming is performed.
#' @param maxit Number of iterations. After zeroing small values, the model is
#'   refitted, and this is repeated until there is nothing to trim or \code{maxit}
#'   iterations are done.
#' @param return_loglik Return the log-likelihood of the trimmed model together with
#'   the model object. The default is \code{FALSE}.
#' @param zerotol Values smaller than this are trimmed to zero.
#' @param verbose Print results of trimming. The default is \code{TRUE}.
#' @param ... Further parameters passed on to \code{\link{fit_model}}.
#'
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_model}} for building and fitting
#' hidden Markov models; and \code{\link{hmm_biofam}} for information on the model used
#' in the example.
#'
#' @examples
#' data("hmm_biofam")
#'
#' # Testing if changing parameter values smaller than 1e-03 to zero
#' # leads to improved log-likelihood.
#' hmm_trim <- trim_model(hmm_biofam, zerotol = 1e-03, maxit = 10)
trim_model <- function(model, maxit = 0, return_loglik=FALSE, zerotol=1e-8, verbose = TRUE, ...){

  ll_original <- logLik(model)
  model_original <- model

  if(inherits(model, "hmm")){

    if(model$n_channels==1){
      if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) ||
          any(model$transition_probs < zerotol & model$transition_probs > 0)
        || any(model$emission_probs < zerotol & model$emission_probs > 0))){
        if(verbose)
          print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }

      model$initial_probs[model$initial_probs < zerotol] <- 0
      model$initial_probs <- model$initial_probs/sum(model$initial_probs)
      model$transition_probs[model$transition_probs < zerotol] <- 0
      model$transition_probs <- model$transition_probs/rowSums(model$transition_probs)
      model$emission_probs[model$emission_probs < zerotol] <- 0
      model$emission_probs <- model$emission_probs/rowSums(model$emission_probs)


      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      if (!isTRUE(all.equal(ll0,ll_original)) && ll0 < ll_original) {
        warning("Trimming resulted model with smaller log-likelihood; returning the original model. ")
        return(model_original)
      }
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_model(model, ...)
          ll <- fit$logLik

          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break

          if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) ||
              any(model$transition_probs < zerotol & model$transition_probs > 0)
            || any(model$emission_probs < zerotol & model$emission_probs > 0)))
            break

          model$initial_probs[model$initial_probs < zerotol] <- 0
          model$initial_probs <- model$initial_probs/sum(model$initial_probs)
          model$transition_probs[model$transition_probs < zerotol] <- 0
          model$transition_probs <- model$transition_probs/rowSums(model$transition_probs)
          model$emission_probs[model$emission_probs < zerotol] <- 0
          model$emission_probs <- model$emission_probs/rowSums(model$emission_probs)

        }
      }

    } else {
      if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) ||
          any(model$transition_probs < zerotol & model$transition_probs > 0)
        || any(sapply(model$emission_probs,function(x) any(x < zerotol & x > 0))))){
        if(verbose)
          print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }
      model$initial_probs[model$initial_probs < zerotol] <- 0
      model$initial_probs <- model$initial_probs/sum(model$initial_probs)
      model$transition_probs[model$transition_probs < zerotol] <- 0
      model$transition_probs <- model$transition_probs/rowSums(model$transition_probs)
      for(i in 1:model$n_channels){
        model$emission_probs[[i]][model$emission_probs[[i]] < zerotol] <- 0
        model$emission_probs[[i]] <- model$emission_probs[[i]]/
          rowSums(model$emission_probs[[i]])

      }


      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      if (!isTRUE(all.equal(ll0,ll_original)) && ll0 < ll_original)  {
        warning("Trimming resulted model with smaller log-likelihood; returning the original model. ")
        return(model_original)
      }
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_model(model, ...)
          ll <- fit$logLik

          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break

          if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) ||
              any(model$transition_probs < zerotol & model$transition_probs > 0)
            || any(sapply(model$emission_probs,function(x) any(x < zerotol & x > 0)))))
            break

          model$initial_probs[model$initial_probs < zerotol] <- 0
          model$initial_probs <- model$initial_probs/sum(model$initial_probs)
          model$transition_probs[model$transition_probs < zerotol] <- 0
          model$transition_probs <- model$transition_probs/rowSums(model$transition_probs)
          for(i in 1:model$n_channels){
            model$emission_probs[[i]][model$emission_probs[[i]] < zerotol] <- 0
            model$emission_probs[[i]] <- model$emission_probs[[i]]/
              rowSums(model$emission_probs[[i]])
          }
        }
      }
    }

  }else if(inherits(model, "mhmm")){
    if(model$n_channels==1){
      if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) ||
          any(unlist(model$transition_probs) < zerotol & unlist(model$transition_probs) > 0)
        || any(unlist(model$emission_probs) < zerotol & unlist(model$emission_probs) > 0))){
        if(verbose)
          print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }

      for(m in 1:model$n_clusters){
        model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
        model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
        model$transition_probs[[m]][model$transition_probs[[m]] < zerotol] <- 0
        model$transition_probs[[m]] <- model$transition_probs[[m]]/rowSums(model$transition_probs[[m]])
        model$emission_probs[[m]][model$emission_probs[[m]] < zerotol] <- 0
        model$emission_probs[[m]] <- model$emission_probs[[m]]/rowSums(model$emission_probs[[m]])
      }

      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      if (!isTRUE(all.equal(ll0,ll_original)) && ll0 < ll_original) {
        warning("Trimming resulted model with smaller log-likelihood.")
        return(model_original)
      }
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_model(model, ...)
          ll <- fit$logLik

          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break

          if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) ||
              any(unlist(model$transition_probs) < zerotol & unlist(model$transition_probs) > 0)
            || any(unlist(model$emission_probs) < zerotol & unlist(model$emission_probs) > 0)))
            break

          for(m in 1:model$n_clusters){
            model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
            model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
            model$transition_probs[[m]][model$transition_probs[[m]] < zerotol] <- 0
            model$transition_probs[[m]] <- model$transition_probs[[m]]/rowSums(model$transition_probs[[m]])
            model$emission_probs[[m]][model$emission_probs[[m]] < zerotol] <- 0
            model$emission_probs[[m]] <- model$emission_probs[[m]]/rowSums(model$emission_probs[[m]])
          }
        }
      }

    } else {
      if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) ||
          any(unlist(model$transition_probs) < zerotol & unlist(model$transition_probs) > 0)
        || any(unlist(model$emission_probs) < zerotol & unlist(model$emission_probs) > 0))){
        if(verbose)
          print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }
      for(m in 1:model$n_clusters){
        model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
        model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
        model$transition_probs[[m]][model$transition_probs[[m]] < zerotol] <- 0
        model$transition_probs[[m]] <- model$transition_probs[[m]]/rowSums(model$transition_probs[[m]])
        for(i in 1:model$n_channels){
          model$emission_probs[[m]][[i]][model$emission_probs[[m]][[i]] < zerotol] <- 0
          model$emission_probs[[m]][[i]] <- model$emission_probs[[m]][[i]]/
            rowSums(model$emission_probs[[m]][[i]])

        }
      }


      if(!is.finite(ll0 <- logLik(model))){
        warning("Trimming resulted in non-finite log-likelihood; returning the original model. Try changing the zerotol parameter.")
        return(model_original)
      }
      if (!isTRUE(all.equal(ll0,ll_original)) && ll0 < ll_original) {
        warning("Trimming resulted model with smaller log-likelihood.")
        return(model_original)
      }
      if(maxit > 0){
        for(ii in 1:maxit){
          fit <- fit_model(model, ...)
          ll <- fit$logLik

          if(ll > ll0){
            model <- fit$model
            ll0 <- ll
          } else break

          if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) ||
              any(unlist(model$transition_probs) < zerotol & unlist(model$transition_probs) > 0)
            || any(unlist(model$emission_probs) < zerotol & unlist(model$emission_probs) > 0)))
            break

          for(m in 1:model$n_clusters){
            model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
            model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
            model$transition_probs[[m]][model$transition_probs[[m]] < zerotol] <- 0
            model$transition_probs[[m]] <- model$transition_probs[[m]]/rowSums(model$transition_probs[[m]])
            for(i in 1:model$n_channels){
              model$emission_probs[[m]][[i]][model$emission_probs[[m]][[i]] < zerotol] <- 0
              model$emission_probs[[m]][[i]] <- model$emission_probs[[m]][[i]]/
                rowSums(model$emission_probs[[m]][[i]])

            }
          }

        }
      }
    }
  }else{
    stop("An object of class hmm or mhmm required.")
  }


  if (verbose) {
    if(maxit > 0)
      print(paste(ii,"iteration(s) used."))

    if(ll0 > ll_original){
      print(paste("Trimming improved log-likelihood, ll_trim-ll_orig =", signif(ll0-ll_original, 3)))
    }
  }


  if(return_loglik){
    list(model=model,loglik=ll0)
  } else model
}
