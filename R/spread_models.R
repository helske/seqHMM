# Transform a mhmm object to separate hmm objects

spread_models <- function(model){

  stnames <- unlist(model$original_state_names)
  rownames(model$transition_probs) <- colnames(model$transition_probs) <- stnames

  if(model$n_channels == 1) {
    rownames(model$emission_probs) <- stnames
  }else{
    for(j in 1:model$n_channels){
      rownames(model$emission_probs[[j]]) <- stnames
    }
    names(model$emission_probs) <- model$channel_names
  }


  transM <- vector("list", model$n_clusters)
  emissM <- vector("list", model$n_clusters)
  init <- vector("list", model$n_clusters)
  k <- 0
  if(model$n_channels == 1){
    for(m in 1:model$n_clusters){
      transM[[m]] <- model$transition_probs[(k+1):(k+model$n_states_in_clusters[m]),
                                            (k+1):(k+model$n_states_in_clusters[m]),drop=FALSE]
        emissM[[m]] <- model$emission_probs[(k+1):(k+model$n_states_in_clusters[m]),,drop=FALSE]
      init[[m]] <- unname(model$initial_probs[(k+1):(k+model$n_states_in_clusters[m])])
      k <- sum(model$n_states_in_clusters[1:m])
    }
  }else{
    for(m in 1:model$n_clusters){
      transM[[m]] <- model$transition_probs[(k+1):(k+model$n_states_in_clusters[m]),
                                            (k+1):(k+model$n_states_in_clusters[m]),drop=FALSE]
      for(j in 1:model$n_channels){
        emissM[[m]][[j]] <- model$emission_probs[[j]][(k+1):(k+model$n_states_in_clusters[m]),,drop=FALSE]
      }
      names(emissM[[m]]) <- model$channel_names
      init[[m]] <- unname(model$initial_probs[(k + 1):(k + model$n_states_in_clusters[m])])
      k <- sum(model$n_states_in_clusters[1:m])
    }
  }

  names(transM) <- names(emissM) <- names(init) <- model$cluster_names


  model$transition_probs <- transM
  model$emission_probs <- emissM
  model$initial_probs <- init
  model$state_names <- model$original_state_names
  model$n_states <- model$n_states_in_clusters
  model$original_state_names <- model$n_states_in_clusters <- NULL
  for (m in 1:model$n_clusters) {
    names(model$initial_probs[[m]]) <- model$state_names[[m]]
  }
  class(model) <- "mhmm"

  model
}
