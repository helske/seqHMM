#' @importFrom Matrix .bdiag
combine_models <- function(model){
  
  n_states_in_clusters <- model$n_states
  n_states <- sum(model$n_states)
  transition_matrix <- as.matrix(.bdiag(model$transition_matrix))
  state_names <- unlist(original_state_names<-model$state_names)
  if (length(unique(state_names))!= length(state_names)){
    state_names <- paste(state_names,rep(1:model$n_clusters,model$n_states),sep="_")
  }
  dimnames(transition_matrix) <- replicate(2, state_names, simplify=FALSE)
  
  if(model$n_channels>1){
    emission_matrix <- lapply(1:model$n_channels, function(i){
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