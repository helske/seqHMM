#' Forward and Backward Probabilities for Hidden Markov Model
#'
#' The \code{forward_backward} function computes scaled forward and backward probabilities of a hidden Markov model.
#'
#' @export 
#' @param model Object of class \code{hmm} or \code{mhmm}.
#' @param forward_only If \code{TRUE}, only forward probabilities are computed. The default is \code{FALSE}.
#' @param log_space Compute forward and backward probabilities in logarithmic scale instead of scaling. 
#'  The default is \code{FALSE}.
#' @param threads Number of threads used in parallel computing. The default is 1.
#' @return List with components 
#'   \item{forward_probs}{If \code{log_space = FALSE}, scaled forward probabilities, i.e. probability of state given 
#'   observations up to that time point. If \code{log_space = TRUE}, 
#'   logarithms of non-scaled forward probabilities. } 
#'   \item{backward_probs}{Scaled backward probabilities (\code{log_space = FALSE}), 
#'   or logarithms of non-scaled backward probabilities(\code{log_space = TRUE}). } 
#'   \item{scaling_factors}{Sum of non-scaled forward probabilities at each time point. 
#'   Only computed if \code{log_space = FALSE}.} 
#'   In case of multiple observations, these are computed independently for each sequence.
#' @examples 
#' # Load a pre-defined MHMM
#' data("mhmm_biofam")
#' 
#' # Compute forward and backward probabilities
#' fb <- forward_backward(mhmm_biofam)
#' 
#' # The most probable hidden state at time t 
#' # given the observations up to time t for the first subject:
#' apply(fb$forward_probs[, , 1], 2, which.max)
#' 
forward_backward <- function(model, forward_only = FALSE, log_space = FALSE, threads = 1){
  
  if(!inherits(model, c("hmm", "mhmm")))
    stop("Argument model must be an object of class 'hmm' or 'mhmm.")
  
  if(inherits(model,"mhmm")){
    mix <- TRUE
    model <- combine_models(model)
  } else mix <- FALSE
  
  if(model$n_channels == 1){
    model$observations <- list(model$observations)
    model$emission_probs <- list(model$emission_probs)
  }
  
  obsArray <- array(0,c(model$n_sequences,model$length_of_sequences,model$n_channels))
  for(i in 1:model$n_channels){
    obsArray[,,i] <- data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$n_symbols[i]] <- model$n_symbols[i]
  }    
  obsArray <- aperm(obsArray)
  emissionArray <- array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
  for(i in 1:model$n_channels)
    emissionArray[,1:model$n_symbols[i],i]<-model$emission_probs[[i]]
  
  if (mix) {
    if (!log_space) {
      out <- forwardbackwardx(model$transition_probs, emissionArray, 
        model$initial_probs, obsArray, model$coefficients,model$X,model$n_states_in_clusters, forward_only, threads)
    } else {
      out <- log_forwardbackwardx(model$transition_probs, emissionArray, 
        model$initial_probs, obsArray, model$coefficients,model$X,model$n_states_in_clusters, forward_only, threads)
    }
  } else{
    if (!log_space) {
      out <- forwardbackward(model$transition_probs, emissionArray, 
        model$initial_probs, obsArray, forward_only, threads)
    } else {
      out <- log_forwardbackward(model$transition_probs, emissionArray, 
        model$initial_probs, obsArray, forward_only, threads)
    }
  }
  
  if (is.null(time_names <- colnames(model$observations[[1]]))) {
    time_names <- 1:model$length_of_sequences
  }
  
  if (is.null(sequence_names <- rownames(model$observations[[1]]))) {
    sequence_names <- 1:model$n_sequences
  } 
  
  dimnames(out$forward_probs) <- list("state" = model$state_names, 
    "time" = time_names, "sequence" = sequence_names)
  if (!log_space) {
    dimnames(out$scaling_factors) <- list("time" = time_names, "sequence" = sequence_names)
  }
  if (!forward_only) {
    dimnames(out$backward_probs) <- list("state" = model$state_names, 
      "time" = time_names, "sequence" = sequence_names)
  }
  
  out
}