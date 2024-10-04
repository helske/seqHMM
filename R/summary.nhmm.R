#' Summary of Non-homogeneous Hidden Markov Models
#'
#' Returns basic information about the model as well as 
#' the estimated regression coefficients.
#'
#' @export
#' @inheritParams get_probs
#' @rdname summary_nhmm
#' @param object Non-homogeneous hidden Markov model of class `nhmm`.
summary.nhmm <- function(object, ...) {
  ll <- logLik(object)
  out <- list(
    logLik = ll, BIC = BIC(ll),
    coefficients = coef(object)
  )
  class(out) <- "summary.nhmm"
  out
}
#'
#' @export
#' @rdname summary_nhmm
#' @param object Non-homogeneous hidden Markov model of class `mnhmm`.
summary.mnhmm <- function(object, ...) {
  cf <- coef(object)
  pr <- exp(object$X_cluster %*% object$gammas$omega)
  prior_cluster_probabilities <- pr / rowSums(pr)
  pcp <- posterior_cluster_probabilities(object)
  mpc <- factor(
    apply(pcp, 1, which.max),
    levels = 1:object$n_clusters, labels = object$cluster_names
  )
  clProbs <- matrix(NA, nrow = object$n_clusters, ncol = object$n_clusters)
  rownames(clProbs) <- colnames(clProbs) <- object$cluster_names
  for (i in 1:object$n_clusters) {
    for (j in 1:object$n_clusters) {
      clProbs[i, j] <- mean(pcp[mpc == object$cluster_names[i], j])
    }
  }
  ll <- logLik(object)
  out <- list(
    logLik = ll, BIC = BIC(ll),
    coefficients = cf,
    most_probable_cluster = mpc,
    prior_cluster_probabilities = prior_cluster_probabilities,
    posterior_cluster_probabilities = pcp,
    classification_table = clProbs
  )
  class(out) <- "summary.mnhmm"
  out
}
#' @export
print.summary.nhmm <- function(x, digits = 3, ...) {
  print(x)
  cat("\nCoefficients:\n")
  print.listof(x$gammas, digits = digits, ...)
  
  cat("Log-likelihood:", x$logLik, "  BIC:", x$BIC, "\n\n")
  invisible(x)
}
#' @export
print.summary.mnhmm <- function(x, digits = 3, ...) {
  print(x)
  cat("\nCoefficients:\n")
  print.listof(x$gammas, digits = digits, ...)
  
  cat("Log-likelihood:", x$logLik, "  BIC:", x$BIC, "\n\n")
  
  cat("Means of prior cluster probabilities :\n")
  print(colMeans(x$prior_cluster_probabilities), digits = digits, ...)
  cat("\n")
  
  tbl <- table(x$most_probable_cluster)
  cl <- matrix(
    c(
      as.character(tbl),
      as.character(round(prop.table(tbl), digits = digits))
    ),
    nrow = 2, byrow = TRUE
  )
  colnames(cl) <- cluster_names
  rownames(cl) <- c("count", "proportion")
  cat("Most probable clusters :\n")
  print.default(cl, quote = FALSE, print.gap = 2, right = TRUE)
  cat("\n")
  
  cat("Classification table :\n")
  cat("Mean cluster probabilities (in columns) by the most probable cluster (rows)\n\n")
  print(x$classification_table, digits = digits, ...)
  invisible(x)
}
