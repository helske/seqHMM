#' Summary method for non-homogeneous hidden Markov models
#'
#' Function \code{summary.nhmm} gives a summary of a non-homogeneous hidden Markov model.
#'
#' @export
#' @method summary nhmm
#' @param object Non-homogeneous hidden Markov model of class \code{nhmm}.
#' @param parameters Whether or not to return transition, emission, and
#' initial probabilities. \code{FALSE} by default.
#' @param conditional_se Return conditional standard errors of coefficients.
#' See \code{\link{vcov.mhmm}} for details. \code{TRUE} by default.
#' @param log_space Make computations using log-space instead of scaling for greater
#' numerical stability at cost of decreased computational performance. Default is \code{FALSE}.
#' @param ... Ignored.
#'
#' @return \describe{
#'    \item{transition_probs}{Transition probabilities. Only returned if \code{parameters = TRUE}.}
#'    \item{emission_probs}{Emission probabilities. Only returned if \code{parameters = TRUE}.}
#'    \item{initial_probs}{Initial state probabilities. Only returned if \code{parameters = TRUE}.}
#'    \item{logLik}{Log-likelihood.}
#'    \item{BIC}{Bayesian information criterion.}
#'    \item{most_probable_cluster}{The most probable cluster according to posterior probabilities.}
#'    \item{coefficients}{Coefficients of covariates.}
#'    \item{vcov}{Variance-covariance matrix of coefficients.}
#'    \item{prior_cluster_probabilities}{Prior cluster probabilities
#'    (mixing proportions) given the covariates.}
#'    \item{posterior_cluster_probabilities}{Posterior cluster membership probabilities.}
#'    \item{classification_table}{Cluster probabilities (columns) by the most probable cluster (rows).}
#'   }
#'
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_model}} for building and
#'   fitting mixture hidden Markov models; and
#'   \code{\link{mhmm_biofam}} for information on the model used in examples.
#'
#' @examples
#' # Loading mixture hidden Markov model (mhmm object)
#' # of the biofam data
#' data("mhmm_biofam")
#'
#' # Model summary
#' summary(mhmm_biofam)
#'
summary.nhmm <- function(object, ...) {


    summary_nhmm <- list(
      transition_probs = object$transition_probs,
      emission_probs = object$emission_probs,
      initial_probs = object$initial_probs,
      logLik = ll, BIC = BIC(ll), most_probable_cluster = most_probable_cluster,
      coefficients = object$coefficients, vcov = vcov(object, conditional_se, log_space = log_space, ...),
      prior_cluster_probabilities = prior_cluster_probabilities,
      posterior_cluster_probabilities = posterior_cluster_probabilities,
      classification_table = clProbs
    )
  
  class(summary_nhmm) <- "summary.nhmm"
  summary_nhmm
}
