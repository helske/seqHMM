#' Hidden Markov model for the biofam data
#' 
#' A hidden Markov model (HMM) fitted for the \code{\link[TraMineR]{biofam}} data.
#' 
#' @format A hidden Markov model of class \code{hmm}; 
#' a left-to-right model with four hidden states.
#'   
#' @seealso Examples of building and fitting HMMs in \code{\link{build_hmm}} and 
#' \code{\link{fit_hmm}}; and \code{\link[TraMineR]{biofam}} for the original data and
#' \code{\link{biofam3c}} for the three-channel version used in this model.
#' 
#' @docType data
#' @keywords datasets
#' @name hmm_biofam
#' @examples
#' data(hmm_biofam)
#' 
#' # Plotting the model
#' plot(hmm_biofam)
#'   
NULL



