#' Hidden Markov model for the biofam data
#'
#' A five-state hidden Markov model (HMM) fitted for the \code{\link[TraMineR]{biofam}} data.
#'
#' @format A hidden Markov model of class \code{hmm};
#' a left-to-right model with four hidden states.
#'
#' @details
#' The model is loaded by calling \code{data(hmm_biofam)}. It was created with the
#' following code:
#' \preformatted{
#' data("biofam3c")
#'
#' # Building sequence objects
#' marr_seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child_seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left_seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' ## Choosing colors
#' attr(marr_seq, "cpal") <- c("violetred2", "darkgoldenrod2", "darkmagenta")
#' attr(child_seq, "cpal") <- c("darkseagreen1", "coral3")
#' attr(left_seq, "cpal") <- c("lightblue", "red3")
#'
#' init <- c(0.9, 0.05, 0.02, 0.02, 0.01)
#'
#' # Starting values for transition matrix
#' trans <- matrix(
#'   c(0.8, 0.10, 0.05, 0.03, 0.02,
#'     0,    0.9, 0.05, 0.03, 0.02,
#'     0,      0,  0.9, 0.07, 0.03,
#'     0,      0,    0,  0.9,  0.1,
#'     0,      0,    0,    0,    1),
#'   nrow = 5, ncol = 5, byrow = TRUE)
#'
#' # Starting values for emission matrices
#' emiss_marr <- matrix(
#'   c(0.9, 0.05, 0.05, # High probability for single
#'     0.9, 0.05, 0.05,
#'     0.05, 0.9, 0.05, # High probability for married
#'     0.05, 0.9, 0.05,
#'     0.3, 0.3, 0.4), # mixed group
#'   nrow = 5, ncol = 3, byrow = TRUE)
#'
#' emiss_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.5, 0.5),
#'   nrow = 5, ncol = 2, byrow = TRUE)
#'
#' emiss_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.5, 0.5),
#'   nrow = 5, ncol = 2, byrow = TRUE)
#'
#' initmod <- build_hmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   initial_probs = init, transition_probs = trans,
#'   emission_probs = list(emiss_marr, emiss_child,
#'     emiss_left),
#'   channel_names = c("Marriage", "Parenthood", "Residence"))
#'
#' fit_biofam <- fit_model(initmod, em = FALSE, local = TRUE)
#' hmm_biofam <- fit_biofam$model
#' }
#'
#' @seealso Examples of building and fitting HMMs in \code{\link{build_hmm}} and
#' \code{\link{fit_model}}; and \code{\link[TraMineR]{biofam}} for the original data and
#' \code{\link{biofam3c}} for the three-channel version used in this model.
#'
#' @docType data
#' @keywords datasets
#' @name hmm_biofam
#' @examples
#'
#' # Plotting the model
#' plot(hmm_biofam)
#'
NULL



