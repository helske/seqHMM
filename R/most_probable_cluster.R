most_probable_cluster <- function(x, type = "viterbi", hp) {
  type <- match.arg(type, c("viterbi", "posterior"))
  if (type == "viterbi") {
    mm <- NULL
    state_names <- unlist(x$state_names)
    clusters <- numeric(x$n_sequences)
    for (i in seq_len(x$n_clusters)) {
      # Find matching cluster names from the first hidden state of each individual
      if (length(unique(state_names)) == length(state_names)) {
        idx <- which(hp[, 1] %in% x$state_names[[i]])
      } else {
        idx <- which(grepl(paste0(x$cluster_names[i], ":"), hp[, 1]))
      }
      clusters[idx] <- x$cluster_names[i]
    }
  } else {
    clusters <- apply(posterior_cluster_probabilities(x), 1, which.max)
  }
  factor(clusters, levels = seq_len(x$n_clusters), labels = x$cluster_names)
}

posterior_cluster_probabilities <- function(x) {
  pp <- posterior_probs(x, as_data_frame = FALSE)
  posterior_cluster_probabilities <- matrix(0, x$n_sequences, x$n_clusters)
  n_states <- rep(x$n_states, length.out = x$n_channels)
  p <- 0
  for (i in seq_len(x$n_clusters)) {
    posterior_cluster_probabilities[, i] <- 
      colSums(pp[(p + 1):(p + n_states[i]), 1, ])
    p <- p + n_states[i]
  }
  dimnames(posterior_cluster_probabilities) <- list(
    sequence = rownames(x$observations), cluster = x$cluster_names
  )
  posterior_cluster_probabilities
}