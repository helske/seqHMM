#' Mixture hidden Markov model for the biofam data
#' 
#' A mixture hidden Markov model (MHMM) fitted for the \code{\link{TraMineR::biofam}} data.
#' 
#' @format A mixture hidden Markov model of class \code{mhmm}: 
#' three clusters with left-to-right models including 4, 4, and 6 hidden states.
#' Two covariates, \code{sex} and \code{cohort}, explaining cluster membership.
#'   
#' @seealso Examples of building and fitting MHMMs in \code{\link{build_mhmm}} and 
#' \code{\link{fit_mhmm}}; and \code{\link{TraMineR::biofam}} for the original data and
#' \code{\link{biofam3c}} for the three-channel version used in this model.
#' 
#' @docType data
#' @keywords datasets
#' @name mhmm_biofam
#' @examples
#' data(mhmm_biofam)
#' 
#' summary(mhmm_biofam)
#' 
#' #' \dontrun{
#' # Plotting the model for each cluster (change with Enter)
#' plot(mhmm_biofam)
#' }
#'   
NULL



