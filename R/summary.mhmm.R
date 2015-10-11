#' Summary method for mixture hidden Markov models
#' 
#' Function \code{summary.mhmm} gives a summary of a mixture hidden Markov model.
#'   
#' 
#' @export
#' @method summary mhmm 
#' @param object Mixture hidden Markov model of class \code{mhmm}.
#' @param parameters Whether or not to print parameters of transition, emission, and 
#' initial probabilities. \code{FALSE} by default.
#' @param conditional_se Return conditional standard errors of coefficients. 
#' See \code{\link{mcov.mhmm}} for details. \code{TRUE} by default.
#' @param ... Further arguments to \code{\link{mcov.mhmm}}.
#' 
#' @details The \code{summary.mhmm} function computes features from a mixture hidden Markov
#' model and stores them as a list. A \code{print} method prints summaries of these:
#' log-likelihood and BIC, coefficients and standard errors of covariates, means of prior 
#' cluster probabilities, and information on most probable clusters.
#' 
#' @return \describe{
#'    \item{logLik}{Log-likelihood}
#'    \item{BIC}{Bayesian information criterion}
#'    \item{coefficients}{Coefficients of covariates}
#'    \item{coef_se}{Standard errors for coefficients}
#'    \item{most_probable_cluster}{The most probable cluster according to posterior probabilities}
#'    \item{prior_cluster_probabilities}{Cluster probabilities (mixing proportions) given by the covariates}
#'    \item{posterior_cluster_probabilities}{Posterior class membership probabilities}
#'    \item{classification_table}{Cluster probabilities (columns) by the most probable cluster (rows)}
#'   }
#'   
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_mhmm}} for building and 
#'   fitting mixture hidden Markov models; and 
#'   \code{\link{mhmm_biofam}} for information on the model used in examples.
#'  
#' @examples 
#' # Loading mixture hidden Markov model (mhmm object)
#' # of the biofam data
#' data(mhmm_biofam)
#'   
#' # Model summary
#' summary(mhmm_biofam)
#'   
#' @seealso \code{\link{fit_mhmm}} for building and fitting mixture hidden Markov models.
#'   

summary.mhmm <- function(object, parameters = FALSE, conditional_se = TRUE, ...){
  
  ll <- logLik(object)
  
  fw <- forward_backward(object, forward_only = TRUE)$forward_probs[,object$length_of_sequences,]
  
  pr <- exp(object$X%*%object$coefficients)
  prior_cluster_probabilities <- pr/rowSums(pr)


  posterior_cluster_probabilities <- prior_cluster_probabilities
  p <- 0
  for(i in 1:object$n_clusters){
    posterior_cluster_probabilities[,i] <- colSums(fw[(p+1):(p+object$n_states[i]), , drop = FALSE])
    p <- p + object$n_states[i]
  }
  most_probable_cluster <- factor(apply(posterior_cluster_probabilities, 1, which.max), 
    levels = 1:object$n_clusters, labels = object$cluster_names)
  

  clProbs <- matrix(NA, nrow = object$n_clusters, ncol = object$n_clusters)
  rownames(clProbs) <- colnames(clProbs) <- object$cluster_names
  for(i in 1:object$n_clusters){
    for(j in 1:object$n_clusters){
      clProbs[i,j] <- mean(posterior_cluster_probabilities[most_probable_cluster == object$cluster_names[i], j])
    }
  }
  
  if(!parameters){
    summary_mhmm <- list(
      logLik = ll, BIC = BIC(ll), most_probable_cluster = most_probable_cluster, 
      coefficients = object$coefficients, vcov = vcov(object, conditional_se, ...),
      prior_cluster_probabilities = prior_cluster_probabilities, 
      posterior_cluster_probabilities = posterior_cluster_probabilities,
      classification_table = clProbs,
      model = object
    )
  }else{
    summary_mhmm <- list(
      transition_matrix = object$transition_matrix,
      emission_matrix = object$emission_matrix,
      initial_probs = object$initial_probs,
      logLik = ll, BIC = BIC(ll), most_probable_cluster = most_probable_cluster, 
      coefficients = object$coefficients, vcov = vcov(object, conditional_se, ...),
      prior_cluster_probabilities = prior_cluster_probabilities, 
      posterior_cluster_probabilities = posterior_cluster_probabilities,
      classification_table = clProbs,
      model = object
    )
  }
  class(summary_mhmm) <- "summary.mhmm"
  summary_mhmm
}