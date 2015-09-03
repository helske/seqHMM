# Transform a mhmm object to separate hmm objects

spread_models <- function(model){
  
  stnames <-unlist(model$original_state_names)
  rownames(model$transition_matrix) <- colnames(model$transition_matrix) <- stnames
  
  if(model$number_of_channels==1){
    rownames(model$emission_matrix) <- stnames
  }else{  
    for(j in 1:model$number_of_channels){
      rownames(model$emission_matrix[[j]]) <- stnames
    }
  }
  
  names(model$emission_matrix) <- model$channel_names
  
  transM <- vector("list", model$number_of_clusters)
  emissM <- vector("list", model$number_of_clusters)
  init <- vector("list", model$number_of_clusters)
  k <- 0
  if(model$number_of_channels==1){
    for(m in 1:model$number_of_clusters){
      transM[[m]] <- model$transition_matrix[(k+1):(k+model$number_of_states_in_clusters[m]),
                                            (k+1):(k+model$number_of_states_in_clusters[m])]
        emissM[[m]] <- model$emission_matrix[(k+1):(k+model$number_of_states_in_clusters[m]),]
      names(emissM[[m]]) <- model$channel_names
      init[[m]] <- unname(model$initial_probs[(k+1):(k+model$number_of_states_in_clusters[m])])
      k <- sum(model$number_of_states_in_clusters[1:m])
    }
  }else{  
    for(m in 1:model$number_of_clusters){
      transM[[m]] <- model$transition_matrix[(k+1):(k+model$number_of_states_in_clusters[m]),
                                            (k+1):(k+model$number_of_states_in_clusters[m])]
      for(j in 1:model$number_of_channels){
        emissM[[m]][[j]] <- model$emission_matrix[[j]][(k+1):(k+model$number_of_states_in_clusters[m]),]
      }
      names(emissM[[m]]) <- model$channel_names
      init[[m]] <- unname(model$initial_probs[(k+1):(k+model$number_of_states_in_clusters[m])])
      k <- sum(model$number_of_states_in_clusters[1:m])
    }
  }
  
  names(transM) <- names(emissM) <- names(init) <- model$cluster_names
  
  model$transition_matrix <- transM
  model$emission_matrix <- emissM
  model$initial_probs <- init
  model$state_names <- model$original_state_names
  model$number_of_states <- model$number_of_states_in_clusters
  model$original_state_names<-model$number_of_states_in_clusters<-NULL
  class(model) <- "mhmm"
  
  model
}