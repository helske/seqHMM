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
      # Give a warning, if no subjects assigned to cluster
      if (length(idx) == 0) {
        mm <- c(mm, i)
      }
    }
    if (length(mm) > 0) {
      warning_(
        "When computing the most probable paths, no subjects were assigned to 
        clusters {x$cluster_names[mm]}")
    }
  } else {
    partial_ll <- logLik(object, partials = TRUE, log_space = log_space)
    fw <- forward_backward(object, forward_only = TRUE, log_space = log_space)$forward_probs[, object$length_of_sequences, ]
    pr <- exp(object$X %*% object$coefficients)
    prior_cluster_probabilities <- pr / rowSums(pr)
    posterior_cluster_probabilities <- array(0, dim = dim(pr))
    if (!log_space) {
      p <- 0
      for (i in 1:object$n_clusters) {
        posterior_cluster_probabilities[, i] <- colSums(fw[(p + 1):(p + object$n_states[i]), , drop = FALSE])
        p <- p + object$n_states[i]
      }
    } else {
      for (j in 1:object$n_sequences) {
        p <- 0
        for (i in 1:object$n_clusters) {
          posterior_cluster_probabilities[j, i] <- exp(logSumExp(fw[(p + 1):(p + object$n_states[i]), j]) - partial_ll[j])
          p <- p + object$n_states[i]
        }
      }
    }
    most_probable_cluster <- factor(apply(posterior_cluster_probabilities, 1, which.max),
                                    levels = 1:object$n_clusters, labels = object$cluster_names
    )
  }
  clusters
}
