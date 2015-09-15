#' Mixture hidden Markov model for the mvad data
#' 
#' A mixture hidden Markov model (MHMM) fitted for the \code{\link[TraMineR]{mvad}} data.
#' 
#' @format A mixture hidden Markov model of class \code{mhmm}: 
#' two clusters including 3 and 4 hidden states.
#' No covariates.
#'   
#' @seealso Examples of building and fitting MHMMs in \code{\link{build_mhmm}} and 
#' \code{\link{fit_mhmm}}; and \code{\link[TraMineR]{mvad}} for more information on the data.
#' 
#' @docType data
#' @keywords datasets
#' @name mhmm_mvad
#' @examples
#' data(mhmm_mvad)
#' 
#' summary(mhmm_mvad)
#' 
#' \dontrun{
#' # Plotting the model for each cluster (change with Enter)
#' plot(mhmm_mvad)
#' }
#'   
NULL



