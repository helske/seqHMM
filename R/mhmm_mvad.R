#' Mixture hidden Markov model for the mvad data
#' 
#' A mixture hidden Markov model (MHMM) fitted for the \code{\link[TraMineR]{mvad}} data.
#' 
#' @format A mixture hidden Markov model of class \code{mhmm}: 
#' two clusters including 3 and 4 hidden states.
#' No covariates.
#' 
#' 
#' @details 
#' The model is loaded by calling \code{data(mhmm_mvad)}. It was created with the 
#' following code:
#' \preformatted{
#' data("mvad", package = "TraMineR")
#' 
#' mvad_alphabet <- 
#'   c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad_labels <- c("employment", "further education", "higher education",
#'   "joblessness", "school", "training")
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet, 
#'   states = mvad_scodes, labels = mvad_labels, xtstep = 6)
#' 
#' attr(mvad_seq, "cpal") <- colorpalette[[6]]
#' 
#' # Starting values for the emission matrices
#' emiss_1 <- matrix(
#'   c(0.01, 0.17, 0.01, 0.01, 0.05, 0.75,
#'     0.95, 0.01, 0.01, 0.01, 0.01, 0.01,
#'     0.01, 0.01, 0.01, 0.95, 0.01, 0.01), 
#'   nrow = 3, ncol = 6, byrow = TRUE)
#' 
#' emiss_2 <- matrix(
#'   c(0.01, 0.01, 0.01, 0.01, 0.95, 0.01,
#'     0.01, 0.84, 0.01, 0.09, 0.01, 0.04,
#'     0.01, 0.01, 0.95, 0.01, 0.01, 0.01,
#'     0.95, 0.01, 0.01, 0.01, 0.01, 0.01), 
#'   nrow = 4, ncol = 6, byrow = TRUE)
#' 
#' # Starting values for the transition matrix
#' 
#' trans_1 <-  matrix(
#'   c(0.92, 0.05, 0.03,
#'     0.01, 0.97, 0.02,
#'     0.02, 0.04, 0.94), 
#'   nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' trans_2 <-  matrix(
#'   c(0.93, 0.02, 0.03, 0.02,
#'     0.01, 0.93, 0.02, 0.04,
#'     0.01, 0.01, 0.96, 0.02,
#'     0.01, 0.02, 0.02, 0.95), 
#'   nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' # Starting values for initial state probabilities
#' initial_probs_1 <- c(0.73, 0.23, 0.04)
#' initial_probs_2 <- c(0.4, 0.05, 0.5, 0.05)
#' 
#' # Building a hidden Markov model with starting values
#' # No covariates
#' init_mhmm_mvad <- build_mhmm(
#'   observations = mvad.seq, 
#'   transition_probs = list(trans_1, trans_2), 
#'   emission_probs = list(emiss_1, emiss_2), 
#'   initial_probs = list(initial_probs_1, initial_probs_2))
#' 
#' # Fit the model
#' set.seed(123)
#' fit_mvad <- fit_model(init_mhmm_mvad, control_em = list(restarts = 10))
#' 
#' # Trim the model
#' mhmm_mvad <- trim_model(fit_mvad$model, zerotol = 1e-04)
#' }
#'   
#' @seealso Examples of building and fitting MHMMs in \code{\link{build_mhmm}} and 
#' \code{\link{fit_model}}; and \code{\link[TraMineR]{mvad}} for more information on the data.
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



