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
#' @seealso \code{\link{build_hmm}} for building Hidden Markov models before 
#'   fitting and \code{\link{fit_hmm}} for fitting Hidden Markov models.
#'   
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <- bf == 4 | bf == 5 | bf == 6
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
#' initial_probs <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = initial_probs
#'   )
#' 
#' # Fitting hidden Markov model
#' HMM <- fit_hmm(bHMM)
#' 
#' # Testing if changing parameter values smaller than 1e-03 to zero 
#' # leads to improved log-likelihood.
#' HMMtrim <- trim_hmm(HMM$model, zerotol=1e-03, maxit=10)
#' 
trim_hmm <- function(model, maxit = 0, return_loglik=FALSE, zerotol=1e-8, ...){
  
  if(inherits(model, "hmm")){
    
    if(model$number_of_channels==1){
      if(!(any(model$initial_probs < zerotol & model$initial_probs > 0) || 
             any(model$transition_matrix < zerotol & model$transition_matrix > 0)
           || any(model$emission_matrix < zerotol & model$emission_matrix > 0))){
        
        print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }
      ll_original <- logLik(model)
      
      model$initial_probs[model$initial_probs < zerotol] <- 0
      model$initial_probs <- model$initial_probs/sum(model$initial_probs)
      model$transition_matrix[model$transition_matrix < zerotol] <- 0
      model$transition_matrix <- model$transition_matrix/rowSums(model$transition_matrix)
      model$emission_matrix[model$emission_matrix < zerotol] <- 0
      model$emission_matrix <- model$emission_matrix/rowSums(model$emission_matrix)
      
      
      if(!is.finite(ll0 <- logLik(model)))
        stop("Initial trimming resulted a non-finite log-likelihood. Try changing the zerotol parameter.")
      
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
      ll_original <- logLik(model)
      
      model$initial_probs[model$initial_probs < zerotol] <- 0
      model$initial_probs <- model$initial_probs/sum(model$initial_probs)
      model$transition_matrix[model$transition_matrix < zerotol] <- 0
      model$transition_matrix <- model$transition_matrix/rowSums(model$transition_matrix)
      for(i in 1:model$number_of_channels){
        model$emission_matrix[[i]][model$emission_matrix[[i]] < zerotol] <- 0
        model$emission_matrix[[i]] <- model$emission_matrix[[i]]/
          rowSums(model$emission_matrix[[i]])
        
      }
      
      
      if(!is.finite(ll0 <- logLik(model)))
        stop("Initial trimming resulted a non-finite log-likelihood. Try changing the zerotol parameter.")
      
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
          for(i in 1:model$number_of_channels){
            model$emission_matrix[[i]][model$emission_matrix[[i]] < zerotol] <- 0
            model$emission_matrix[[i]] <- model$emission_matrix[[i]]/
              rowSums(model$emission_matrix[[i]])          
          }       
        }
      }
    }
    
  }else if(inherits(model, "mhmm")){
    if(model$number_of_channels==1){
      if(!(any(unlist(model$initial_probs) < zerotol & unlist(model$initial_probs) > 0) || 
             any(unlist(model$transition_matrix) < zerotol & unlist(model$transition_matrix) > 0)
           || any(unlist(model$emission_matrix) < zerotol & unlist(model$emission_matrix) > 0))){
        
        print("Nothing to trim.")
        if(return_loglik){
          return(list(model=model,loglik=logLik(model)))
        } else return(model)
      }
      ll_original <- logLik(model)
      
      for(m in 1:model$number_of_clusters){
        model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
        model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
        model$transition_matrix[[m]][model$transition_matrix[[m]] < zerotol] <- 0
        model$transition_matrix[[m]] <- model$transition_matrix[[m]]/rowSums(model$transition_matrix[[m]])
        model$emission_matrix[[m]][model$emission_matrix[[m]] < zerotol] <- 0
        model$emission_matrix[[m]] <- model$emission_matrix[[m]]/rowSums(model$emission_matrix[[m]])
      }
      
      if(!is.finite(ll0 <- logLik(model)))
        stop("Initial trimming resulted a non-finite log-likelihood. Try changing the zerotol parameter.")
      
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
          
          for(m in 1:model$number_of_clusters){
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
      ll_original <- logLik(model)
      
      for(m in 1:model$number_of_clusters){
        model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
        model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
        model$transition_matrix[[m]][model$transition_matrix[[m]] < zerotol] <- 0
        model$transition_matrix[[m]] <- model$transition_matrix[[m]]/rowSums(model$transition_matrix[[m]])
        for(i in 1:model$number_of_channels){
          model$emission_matrix[[m]][[i]][model$emission_matrix[[m]][[i]] < zerotol] <- 0
          model$emission_matrix[[m]][[i]] <- model$emission_matrix[[m]][[i]]/
            rowSums(model$emission_matrix[[m]][[i]])
          
        }
      }
      
      
      if(!is.finite(ll0 <- logLik(model)))
        stop("Initial trimming resulted a non-finite log-likelihood. Try changing the zerotol parameter.")
      
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
          
          for(m in 1:model$number_of_clusters){
            model$initial_probs[[m]][model$initial_probs[[m]] < zerotol] <- 0
            model$initial_probs[[m]] <- model$initial_probs[[m]]/sum(model$initial_probs[[m]])
            model$transition_matrix[[m]][model$transition_matrix[[m]] < zerotol] <- 0
            model$transition_matrix[[m]] <- model$transition_matrix[[m]]/rowSums(model$transition_matrix[[m]])
            for(i in 1:model$number_of_channels){
              model$emission_matrix[[m]][[i]][model$emission_matrix[[m]][[i]] < zerotol] <- 0
              model$emission_matrix[[m]][[i]] <- model$emission_matrix[[m]][[i]]/
                rowSums(model$emission_matrix[[m]][[i]])
              
            }
          }
          
        }
      }
    }
  }else{
    stop("An object of class hmm or mixhmm required.")
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