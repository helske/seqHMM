#' Summary method for mixture hidden Markov models
#'
#' Function `summary.mhmm` gives a summary of a mixture hidden Markov model.
#'
#' @export
#' @param object Mixture hidden Markov model of class `mhmm`.
#' @param parameters Whether or not to return transition, emission, and
#' initial probabilities. `FALSE` by default.
#' @param conditional_se Return conditional standard errors of coefficients.
#' See [vcov.mhmm()] for details. `TRUE` by default.
#' @param ... Further arguments to [vcov.mhmm()].
#'
#' @details The `summary.mhmm` function computes features from a mixture hidden Markov
#' model and stores them as a list. A `print` method prints summaries of these:
#' log-likelihood and BIC, coefficients and standard errors of covariates, means of prior
#' cluster probabilities, and information on most probable clusters.
#'
#' @return
#' * transition_probs\cr Transition probabilities. Only returned if `parameters = TRUE`.
#' * emission_probs\cr Emission probabilities. Only returned if `parameters = TRUE`.
#' * initial_probs\cr Initial state probabilities. Only returned if `parameters = TRUE`.
#' * logLik\cr Log-likelihood.
#' * BIC\cr Bayesian information criterion.
#' * most_probable_cluster\cr The most probable cluster according to posterior probabilities.
#' * coefficients\cr Coefficients of covariates.
#' * vcov\cr Variance-covariance matrix of coefficients.
#' * prior_cluster_probabilities\cr Prior cluster probabilities
#'    (mixing proportions) given the covariates.
#' * posterior_cluster_probabilities\cr Posterior cluster membership probabilities.
#' * classification_table\cr Cluster probabilities (columns) by the most probable cluster (rows).
#'
#'
#' @seealso [build_mhmm()] and [fit_model()] for building and
#'   fitting mixture hidden Markov models; and
#'   [mhmm_biofam()] for information on the model used in examples.
#'
#' @examples
#' # Loading mixture hidden Markov model (mhmm object)
#' # of the biofam data
#' data("mhmm_biofam")
#'
#' # Model summary
#' summary(mhmm_biofam)
#'
summary.mhmm <- function(
    object, parameters = FALSE, conditional_se = TRUE, ...) {
  # avoid CRAN check warning due to NSE
  probability <- id <- cluster <- NULL
  pcp <- posterior_cluster_probabilities(object)
  pr <- exp(object$X %*% object$coefficients)
  prior_cluster_probabilities <- data.table(
    pcp[, 1:2], 
    probability = c(t(pr / rowSums(pr)))
  )
  mpc <- pcp[, .SD[which.max(probability)], by = id]$cluster
  clProbs <- matrix(NA, nrow = object$n_clusters, ncol = object$n_clusters)
  rownames(clProbs) <- colnames(clProbs) <- object$cluster_names
  for (i in object$cluster_names) {
    for (j in object$cluster_names) {
      clProbs[j, i] <- mean(pcp[cluster == i][mpc == j, probability])
    }
  }
  ll <- logLik(object)
  if (!parameters) {
    summary_mhmm <- list(
      logLik = ll, BIC = stats::BIC(ll), most_probable_cluster = mpc,
      coefficients = object$coefficients, 
      vcov = vcov(object, conditional_se, ...),
      prior_cluster_probabilities = prior_cluster_probabilities,
      posterior_cluster_probabilities = pcp,
      classification_table = clProbs
    )
  } else {
    summary_mhmm <- list(
      transition_probs = object$transition_probs,
      emission_probs = object$emission_probs,
      initial_probs = object$initial_probs,
      logLik = ll, BIC = stats::BIC(ll), most_probable_cluster = mpc,
      coefficients = object$coefficients, 
      vcov = vcov(object, conditional_se, ...),
      prior_cluster_probabilities = prior_cluster_probabilities,
      posterior_cluster_probabilities = pcp,
      classification_table = clProbs
    )
  }
  class(summary_mhmm) <- "summary_mhmm"
  summary_mhmm
}
#' @export
#' @rdname print
print.summary_mhmm <- function(x, digits = 3, ...) {
  # avoid CRAN check warnings due to NSE
  probability <- cluster <- NULL
  
  if (exists("transition_probs", x)) {
    cat("Initial probabilities :\n")
    print.listof(x$initial_probs, digits = digits, ...)
    cat("Transition probabilities :\n")
    print.listof(x$transition_probs, digits = digits, ...)
    cat("Emission probabilities :\n")
    if (!is.list(x$emission_probs[[1]])) {
      print.listof(x$emission_probs, digits = digits, ...)
    } else {
      for (i in 1:length(x$emission_probs)) {
        cat(names(x$emission_probs)[i], ":\n\n")
        print.listof(x$emission_probs[[i]], digits = digits, ...)
      }
    }
  }
  coef_se <- matrix(sqrt(diag(x$vcov)), nrow(x$coefficients))
  coefs <- replicate((ncol(x$coefficients) - 1),
                     matrix(NA, nrow = nrow(x$coefficients), ncol = 2),
                     simplify = FALSE
  )
  for (i in 1:length(coefs)) {
    coefs[[i]][, 1] <- x$coefficients[, i + 1]
    coefs[[i]][, 2] <- coef_se[, i]
    rownames(coefs[[i]]) <- rownames(x$coefficients)
    colnames(coefs[[i]]) <- c("Estimate", "Std. error")
  }
  cluster_names <- colnames(x$coefficients)
  names(coefs) <- cluster_names[-1]
  cat("Covariate effects :\n")
  cat(cluster_names[1], "is the reference.\n\n")
  print.listof(coefs, print.gap = 2, digits = digits, quote = FALSE, ...)
  
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

