#' Summary method for mixture non-homogenous hidden Markov models
#'
#' @export
#' @param object Non-homogeneous hidden Markov model of class `mnhmm`.
#' @param ... Ignored
summary.mnhmm <- function(object, ...) {
  cf <- coef(object)
  pr <- exp(object$X_omega %*% object$gammas$gamma_omega)
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
    model = object,
    logLik = ll, BIC = stats::BIC(ll),
    coefficients = cf,
    most_probable_cluster = mpc,
    prior_cluster_probabilities = prior_cluster_probabilities,
    posterior_cluster_probabilities = pcp,
    classification_table = clProbs
  )
  class(out) <- "summary_mnhmm"
  out
}
#' @export
print.summary_mnhmm <- function(x, digits = 3, ...) {
  # to avoid NSE warnings
  probability <- cluster <- NULL
  print(x$model)
  
  cat("Log-likelihood:", x$logLik, "  BIC:", x$BIC, "\n\n")
  
  cat("Means of prior cluster probabilities :\n")
  p <- x$prior_cluster_probabilities[, list(mean = mean(probability)), by = cluster]
  p <- stats::setNames(p$mean, p$cluster)
  print(p, digits = digits, ...)
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
