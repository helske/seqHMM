#' Most Probable Paths of Hidden States
#' 
#' Function \code{hidden_paths} computes the most probable path of
#' hidden states of a (mixture) hidden Markov model given the observed sequences.
#' 
#' @export
#' @param model Hidden Markov model of class \code{hmm} or
#'  mixture HMM of class \code{mhmm}.
#' 
#' @return Most probable paths of hidden states as an \code{stslist} object
#' (see \code{\link{seqdef}}). The log-probability included as an attribute.
#'   
#' @examples 
#' # Loading a hidden Markov model of the biofam data (hmm object)
#' data(hmm_biofam)
#' 
#' # Computing the most probable paths of hidden states
#' mpp <- hidden_paths(hmm_biofam)
#'   
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting Hidden Markov models; \code{\link{build_mhmm}} and 
#'   \code{\link{fit_mhmm}} for building and fitting mixture hidden Markov models; 
#'   \code{\link{hmm_biofam}} for information on the model used in the example;
#'   and \code{\link{seqIplot}}, \code{\link{ssplot}}, or \code{\link{mssplot}}
#'   for plotting the most probable paths of hidden states.
#'   

hidden_paths <- function(model){
  
  ll <- logLik(model, partials = TRUE)
  if(inherits(model,"mhmm")){
    model <- combine_models(model)
    mix <- TRUE
  } else mix <- FALSE
  
  
  if(model$n_channels == 1){
    model$observations <- list(model$observations)
    model$emission_matrix <- list(model$emission_matrix)
  }
  
  
  model$initial_probs <- log(model$initial_probs)
  model$transition_matrix <- log(model$transition_matrix)
  
  obsArray <- array(0,c(model$n_sequences,model$length_of_sequences,model$n_channels))
  for(i in 1:model$n_channels){
    obsArray[,,i] <- data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$n_symbols[i]] <- model$n_symbols[i]
  }       
  storage.mode(obsArray) <- "integer"
  
  emissionArray <- array(0,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
  for(i in 1:model$n_channels)
    emissionArray[,1:model$n_symbols[i],i] <- log(model$emission_matrix[[i]])
  
  if(mix){
    out <- viterbix(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray, model$coefficients, 
      model$X, model$n_states_in_clusters)
  } else{
    out <- viterbi(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray)
  }
  
  
  if(model$n_sequences==1){
    mpp <- t(model$state_names[out$q+1])
  }else{
    mpp <- apply(out$q+1,2,function(x) model$state_names[x])
  }
  mpp <- suppressWarnings(
    suppressMessages(
      seqdef(
        mpp,alphabet=model$state_names, id=rownames(model$obs[[1]]),
        start=attr(model$obs[[1]],"start"), xtstep=attr(model$obs[[1]],"xtstep")
      )
    )
  )
  
  if (sum(model$n_states) <= 200) {
    attr(mpp, "cpal") <- seqHMM::colorpalette[[sum(model$n_states)]]
  } else {
    cp <- NULL
    k <- 200
    p <- 0
    while(sum(model$n_states) - p > 0){
      cp <- c(cp, seqHMM::colorpalette[[k]])
      p <- p + k
      k <- k - 1
    }
    attr(mpp, "cpal") <- cp[1:sum(model$n_states)]
  }
  
  attr(mpp, "log_prob") <- out$logp
  
  mpp
}
