
combine_models <- function(model){
  
  n_states_in_clusters <- model$n_states
  n_states <- sum(model$n_states)
  transition_probs <- as.matrix(.bdiag(model$transition_probs))
  state_names <- unlist(original_state_names<-model$state_names)
  if (length(unique(state_names))!= length(state_names)){
    state_names <- paste(rep(model$cluster_names, model$n_states), state_names, sep=":")
  }
  dimnames(transition_probs) <- replicate(2, state_names, simplify=FALSE)
  
  if(model$n_channels>1){
    emission_probs <- lapply(1:model$n_channels, function(i){
      x <- do.call("rbind",sapply(model$emission_probs,"[",i))
      rownames(x) <- state_names
      x
    })
    names(emission_probs) <- model$channel_names
  } else {
    emission_probs <- do.call("rbind",model$emission_probs)
    rownames(emission_probs) <- state_names
  }
  
  model <- list(observations = model$observations,
                transition_probs = transition_probs,
                emission_probs=emission_probs,
                initial_probs = unlist(model$initial_probs),
    coefficients = model$coefficients, X = model$X,
                cluster_names=model$cluster_names,
                state_names=state_names,
                symbol_names=model$symbol_names,
                channel_names=model$channel_names,
                length_of_sequences=model$length_of_sequences,
                n_sequences=model$n_sequences,
                n_symbols=model$n_symbols,
                n_states=n_states,
                n_channels=model$n_channels,
                n_covariates=model$n_covariates,
                n_clusters=model$n_clusters,
                n_states_in_clusters=n_states_in_clusters,
                original_state_names = original_state_names)
  class(model)<-"combined_mhmm"
  model
}