combine_models <- function(model){
  
  number_of_states_in_clusters <- model$number_of_states
  number_of_states <- sum(model$number_of_states)
  transition_matrix <- as.matrix(.bdiag(model$transition_matrix))
  state_names <- unlist(original_state_names<-model$state_names)
  if (length(unique(state_names))!= length(state_names)){
    state_names <- paste(state_names,rep(1:model$number_of_clusters,model$number_of_states),sep="_")
  }
  dimnames(transition_matrix) <- replicate(2, state_names, simplify=FALSE)
  
  if(model$number_of_channels>1){
    emission_matrix <- lapply(1:model$number_of_channels, function(i){
      x <- do.call("rbind",sapply(model$emission_matrix,"[",i))
      rownames(x) <- state_names
      x
    })
    names(emission_matrix) <- model$channel_names
  } else {
    emission_matrix <- do.call("rbind",model$emission_matrix)
    rownames(emission_matrix) <- state_names
  }
  
  model <- list(observations = model$observations,
                transition_matrix = transition_matrix,
                emission_matrix=emission_matrix,
                initial_probs = unlist(model$initial_probs),
                beta=model$beta, X=model$X,
                cluster_names=model$cluster_names,
                state_names=state_names,
                symbol_names=model$symbol_names,
                channel_names=model$channel_names,
                length_of_sequences=model$length_of_sequences,
                number_of_sequences=model$number_of_sequences,
                number_of_symbols=model$number_of_symbols,
                number_of_states=number_of_states,
                number_of_channels=model$number_of_channels,
                number_of_covariates=model$number_of_covariates,
                number_of_clusters=model$number_of_clusters,
                number_of_states_in_clusters=number_of_states_in_clusters,
                original_state_names = original_state_names)
  class(model)<-"combined_mhmm"
  model
}