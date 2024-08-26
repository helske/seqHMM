#' Summary of Non-homogeneous Hidden Markov Models
#'
#' Function `summary.nhmm` returns basic information about the model as well as 
#' the estimated regression coefficients.
#'
#' @export
#' @param object Non-homogeneous hidden Markov model of class `nhmm`.
summary.nhmm <- function(object, nsim = 0, probs = c(0.025, 0.975), ...) {
  print(object)
  out <- coef(object, nsim, probs)
  cat("\nCoefficients:\n")
  print.listof(out, digits = digits, ...)
}
#' Summary method for non-homogeneous hidden Markov models
#'
#' Function `summary.nhmm` returns basic information about the model as well as 
#' the estimated regression coefficients.
#'
#' @export
#' @param object Non-homogeneous hidden Markov model of class `mnhmm`.
summary.mnhmm <- function(object, nsim = 0, probs = c(0.025, 0.975), ...) {
  print(object)
  out <- coef(object, nsim, probs)
  cat("\nCoefficients:\n")
  print.listof(out, digits = digits, ...)
  
  pr <- exp(object$X_cluster %*% out)
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
  clProbs <- matrix(NA, nrow = object$n_clusters, ncol = object$n_clusters)
  rownames(clProbs) <- colnames(clProbs) <- object$cluster_names
  for (i in 1:object$n_clusters) {
    for (j in 1:object$n_clusters) {
      clProbs[i, j] <- mean(posterior_cluster_probabilities[most_probable_cluster == object$cluster_names[i], j])
    }
  }
}