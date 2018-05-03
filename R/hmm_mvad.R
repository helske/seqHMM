#' Hidden Markov model for the mvad data
#'
#' A hidden Markov model (MMM) fitted for the \code{\link[TraMineR]{mvad}} data.
#'
#' @format A hidden Markov model of class \code{hmm};
#' unrestricted model with six hidden states.
#'
#' @details
#' Model was created with the following code:
#' \preformatted{
#'
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
#' # Starting values for the emission matrix
#' emiss <- matrix(
#'   c(0.05, 0.05, 0.05, 0.05, 0.75, 0.05, # SC
#'     0.05, 0.75, 0.05, 0.05, 0.05, 0.05, # FE
#'     0.05, 0.05, 0.05, 0.4,  0.05, 0.4,  # JL, TR
#'     0.05, 0.05, 0.75, 0.05, 0.05, 0.05, # HE
#'     0.75, 0.05, 0.05, 0.05, 0.05, 0.05),# EM
#'   nrow = 5, ncol = 6, byrow = TRUE)
#'
#' # Starting values for the transition matrix
#' trans <- matrix(0.025, 5, 5)
#' diag(trans) <- 0.9
#'
#' # Starting values for initial state probabilities
#' initial_probs <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#'
#' # Building a hidden Markov model
#' init_hmm_mvad <- build_hmm(observations = mvad_seq,
#'   transition_probs = trans, emission_probs = emiss,
#'   initial_probs = initial_probs)
#'
#' set.seed(21)
#' fit_hmm_mvad <- fit_model(init_hmm_mvad, control_em = list(restart = list(times = 100)))
#' hmm_mvad <- fit_hmm_mvad$model
#' }
#'
#' @seealso Examples of building and fitting HMMs in \code{\link{build_hmm}} and
#' \code{\link{fit_model}}; and \code{\link[TraMineR]{mvad}} for more information on the data.
#'
#' @docType data
#' @keywords datasets
#' @name hmm_mvad
#' @examples
#' data("hmm_mvad")
#'
#' # Plotting the model
#' plot(hmm_mvad)
#'
NULL



