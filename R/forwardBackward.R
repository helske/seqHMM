#' Forward Probabilities for Hidden Markov Model
#'
#' Function \code{forward_probs} computes the forward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{hmm}.
#' @return Forward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
forward_probs<-function(model, test = FALSE){
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
    out<-forwardx(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray, model$coefficients,model$X,model$n_states_in_clusters)
  } else{
    if(test){
      out<-forward2(model$transition_matrix, emissionArray, 
        model$initial_probs, obsArray)
    } else {
    out<-forward(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray)
    }
  } 
  
  dimnames(out)<-list("state" = model$state_names,"time" = 1:model$length_of_sequences)
  out
}


#' Backward Probabilities for Hidden Markov Model
#'
#' Function \code{backward_probs} computes the backward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{hmm}.
#' @return Backward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
backward_probs<-function(model){
  if(inherits(model,"mhmm"))
    model <- combine_models(model)
  
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
  
  out<-backward(model$transition_matrix, emissionArray, obsArray)
  
  dimnames(out)<-list("state" = model$state_names,"time" = 1:model$length_of_sequences)
  out
  
}
