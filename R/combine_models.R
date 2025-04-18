#' Combine mixture HMM for a single HMM
#'
#' This function is used internally within various functions dealing with `mhmm` 
#' objects (e.g., [fit_model()]).
#' 
#' @noRd
.combine_models <- function(model) {
  n_states_in_clusters <- model$n_states
  n_states <- sum(model$n_states)
  transition_probs <- as.matrix(Matrix::.bdiag(model$transition_probs))
  state_names <- unlist(original_state_names <- model$state_names)
  if (n_unique(state_names) != length(state_names)) {
    state_names <- paste0(rep(model$cluster_names, model$n_states), 
                          ": ", state_names)
  }
  dimnames(transition_probs) <- replicate(2, state_names, simplify = FALSE)

  if (model$n_channels > 1) {
    emission_probs <- lapply(seq_len(model$n_channels), \(i) {
      x <- do.call(rbind, sapply(model$emission_probs, "[", i))
      rownames(x) <- state_names
      x
    })
    names(emission_probs) <- model$channel_names
  } else {
    emission_probs <- do.call(rbind, model$emission_probs)
    rownames(emission_probs) <- state_names
  }
  model$initial_probs <- unlist(model$initial_probs)
  model$transition_probs <- transition_probs
  model$emission_probs <- emission_probs
  model$state_names <- state_names
  model$n_states <- n_states
  model$n_states_in_clusters <- n_states_in_clusters
  model$original_state_names = original_state_names
  class(model) <- "combined_mhmm"
  model
}
