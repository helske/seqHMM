#' Forward Probabilities for Hidden Markov Model
#'
#' Function \code{forward_probs} computes the forward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{hmm}.
#' @return Forward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
forward_backward <- function(model, forward_only = FALSE){
  if(inherits(model,"mhmm")){
    mix <- TRUE
    model <- combine_models(model)
  } else mix <- FALSE
  
  if(model$n_channels == 1){
    model$observations <- list(model$observations)
    model$emission_matrix <- list(model$emission_matrix)
  }
  
  obsArray<-array(0,c(model$n_sequences,model$length_of_sequences,model$n_channels))
  for(i in 1:model$n_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$n_symbols[i]]<-model$n_symbols[i]
  }    
  storage.mode(obsArray)<-"integer"
  emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
  for(i in 1:model$n_channels)
    emissionArray[,1:model$n_symbols[i],i]<-model$emission_matrix[[i]]
  
  if(mix){
    out<-forwardbackwardx(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray, model$coefficients,model$X,model$n_states_in_clusters, forward_only)
  } else{
    out <- forwardbackward(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray, forward_only)
  }


dimnames(out$forward_probs)<-list("state" = model$state_names,"time" = 1:model$length_of_sequences)
if (!forward_only) {
  dimnames(out$forward_probs)<-list("state" = model$state_names,"time" = 1:model$length_of_sequences)
}

out
}