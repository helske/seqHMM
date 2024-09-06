#' Extract Most Probable Cluster for Each Sequence
#' 
#' @param x An object of class `mhmm` or `mnhmm`.
#' @param type A character string specifying the method to use. Either
#' `"viterbi"` (default) or `"posterior"`. Former uses the most probable hidden 
#' path to determine the cluster membership for each sequence, while the latter 
#' finds the cluster which has the largest sum of posterior probabilities of 
#' states of that cluster. 
#' @param hp An output from [hidden_paths()] function. Only used in case of 
#' `type = "viterbi"`. If missing, hidden paths will be computed using `x`.
#' @return A vector containing the most probable cluster for each sequence.
#' @export
most_probable_cluster <- function(x, type = "viterbi", hp = NULL) {
  stopifnot_(
    inherits(x, "mhmm") || inherits(x, "mnhmm"),
    "Argument {.arg x} must be a {.cls mhmm} or {.cls mnhmm} object."
  )
  type <- match.arg(type, c("viterbi", "posterior"))
  if (type == "viterbi") {
    if (is.null(hp)) hp <- hidden_paths(x)
    mm <- NULL
    state_names <- unlist(x$state_names)
    clusters <- numeric(x$n_sequences)
    for (i in seq_len(x$n_clusters)) {
      # Find matching cluster names from the first hidden state of each individual
      if (length(unique(state_names)) == length(state_names)) {
        idx <- which(hp[, 1] %in% x$state_names[[i]])
      } else {
        idx <- which(grepl(paste0(x$cluster_names[i], ": "), hp[, 1]))
      }
      clusters[idx] <- i
    }
  } else {
    clusters <- apply(posterior_cluster_probabilities(x), 1, which.max)
  }
  factor(clusters, levels = seq_len(x$n_clusters), labels = x$cluster_names)
}
#' Extract Posterior Cluster Probabilities
#' 
#' @param x An object of class `mhmm` or `mnhmm`.
#' @return matrix of posterior cluster probabilities for each sequence and 
#' cluster.
#' @export
posterior_cluster_probabilities <- function(x) {
  stopifnot_(
    inherits(x, "mhmm") || inherits(x, "mnhmm"),
    "Argument {.arg x} must be a {.cls mhmm} or {.cls mnhmm} object."
  )
  pp <- posterior_probs(x, as_data_frame = FALSE)
  posterior_cluster_probabilities <- matrix(0, x$n_sequences, x$n_clusters)
  n_states <- rep(x$n_states, length.out = x$n_clusters)
  p <- 0
  for (i in seq_len(x$n_clusters)) {
    posterior_cluster_probabilities[, i] <- 
      colSums(array(
        pp[(p + 1):(p + n_states[i]), 1, ], 
        dim = c(n_states[i], x$n_sequences)
      ))
    p <- p + n_states[i]
  }
  dimnames(posterior_cluster_probabilities) <- list(
    sequence = rownames(x$observations), cluster = x$cluster_names
  )
  posterior_cluster_probabilities
}