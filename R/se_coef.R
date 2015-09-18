#' Standard Errors for Regression Coefficients of Mixture Hidden Markov Model
#'
#' @param model Object of class \code{mhmm}.
#'
#' @return Matrix containing the standard errors for coefficients.
#' @export
#'
se_coef <- function(model){
  matrix(c(rep(0,model$n_covariates),
    sqrt(diag(varcoef(model$coefficients, model$X, model$n_states)))),
    nrow = model$n_covariates, ncol = model$n_clusters)
}