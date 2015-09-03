#' Forward Probabilities for Hidden Markov Model
#'
#' Function \code{forwardProbs} computes the forward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{hmm}.
#' @return Forward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
forward_probs<-function(model){
  if(inherits(model,"mhmm")){
    mix <- TRUE
    model <- combine_models(model)
  } else mix <- FALSE
  
  if(model$number_of_channels == 1){
    model$observations <- list(model$observations)
    model$emission_matrix <- list(model$emission_matrix)
  }
  
  obsArray<-array(0,c(model$number_of_sequences,model$length_of_sequences,model$number_of_channels))
  for(i in 1:model$number_of_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$number_of_symbols[i]]<-model$number_of_symbols[i]
  }    
  storage.mode(obsArray)<-"integer"
  emissionArray<-array(1,c(model$number_of_states,max(model$number_of_symbols)+1,model$number_of_channels))
  for(i in 1:model$number_of_channels)
    emissionArray[,1:model$number_of_symbols[i],i]<-model$emission_matrix[[i]]
  if(mix){
    out<-forwardx(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray, model$beta,model$X,model$number_of_states_in_clusters)
  } else{
    out<-forward(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray)
  } 
  
  dimnames(out)<-list("state" = model$state_names,"time" = 1:model$length_of_sequences)
  out
}


#' Backward Probabilities for Hidden Markov Model
#'
#' Function \code{backwardProbs} computes the backward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{hmm}.
#' @return Backward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
backward_probs<-function(model){
  if(inherits(model,"mhmm"))
    model <- combine_models(model)
  
  if(model$number_of_channels == 1){
    model$observations <- list(model$observations)
    model$emission_matrix <- list(model$emission_matrix)
  }
  
  obsArray<-array(0,c(model$number_of_sequences,model$length_of_sequences,model$number_of_channels))
  for(i in 1:model$number_of_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$number_of_symbols[i]]<-model$number_of_symbols[i]
  }    
  storage.mode(obsArray)<-"integer"
  emissionArray<-array(1,c(model$number_of_states,max(model$number_of_symbols)+1,model$number_of_channels))
  for(i in 1:model$number_of_channels)
    emissionArray[,1:model$number_of_symbols[i],i]<-model$emission_matrix[[i]]
  
  out<-backward(model$transition_matrix, emissionArray, obsArray)
  
  dimnames(out)<-list("state" = model$state_names,"time" = 1:model$length_of_sequences)
  out
  
}
