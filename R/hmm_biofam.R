#' Hidden Markov model for the biofam data
#' 
#' A hidden Markov model (HMM) fitted for the \code{\link[TraMineR]{biofam}} data.
#' 
#' @format A hidden Markov model of class \code{hmm}; 
#' a left-to-right model with four hidden states.
#' 
#' @details 
#' The model is loaded by calling \code{data(hmm_biofam)}. It was created with the 
#' following code:
#' \preformatted{
#' data(biofam3c)
#' 
#' # Building sequence objects
#' child.seq <- seqdef(biofam3c$children, start = 15)
#' marr.seq <- seqdef(biofam3c$married, start = 15)
#' left.seq <- seqdef(biofam3c$left, start = 15)
#' 
#' ## Choosing colors
#' attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#' 
#' # Starting values for emission matrices
#' emiss_marr <- matrix(NA, nrow=4, ncol=3)
#' emiss_marr[1,] <- seqstatf(marr.seq[, 1:4])[, 2] + 0.1
#' emiss_marr[2,] <- seqstatf(marr.seq[, 5:8])[, 2] + 0.1
#' emiss_marr[3,] <- seqstatf(marr.seq[, 9:12])[, 2] + 0.1
#' emiss_marr[4,] <- seqstatf(marr.seq[, 13:16])[, 2] + 0.1
#' emiss_marr <- emiss_marr / rowSums(emiss_marr)
#' 
#' emiss_child <- matrix(NA, nrow=4, ncol=2)
#' emiss_child[1,] <- seqstatf(child.seq[, 1:4])[, 2] + 0.1
#' emiss_child[2,] <- seqstatf(child.seq[, 5:8])[, 2] + 0.1
#' emiss_child[3,] <- seqstatf(child.seq[, 9:12])[, 2] + 0.1
#' emiss_child[4,] <- seqstatf(child.seq[, 13:16])[, 2] + 0.1
#' emiss_child <- emiss_child / rowSums(emiss_child)
#' 
#' emiss_left <- matrix(NA, nrow=4, ncol=2)
#' emiss_left[1,] <- seqstatf(left.seq[, 1:4])[, 2] + 0.1
#' emiss_left[2,] <- seqstatf(left.seq[, 5:8])[, 2] + 0.1
#' emiss_left[3,] <- seqstatf(left.seq[, 9:12])[, 2] + 0.1
#' emiss_left[4,] <- seqstatf(left.seq[, 13:16])[, 2] + 0.1
#' emiss_left <- emiss_left / rowSums(emiss_left)
#' 
#' # Initial values for transition matrix
#' trans <- matrix(
#'   c(0.9, 0.06, 0.03, 0.01,
#'       0,    0.9, 0.07, 0.03,
#'       0,      0,  0.9,  0.1,
#'       0,      0,    0,    1), 
#'   nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities
#' initial_probs <- c(0.9, 0.07, 0.02, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' init_hmm_biofam <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_matrix = trans,
#'   emission_matrix = list(emiss_child, emiss_marr, emiss_left),
#'   initial_probs = initial_probs,
#'   channel_names = c("Parenthood", "Marriage", "Residence"),
#'   state_names = paste("State", 1:4))
#' 
#' fit_biofam <- fit_hmm(
#'    init_hmm_biofam, control_global = list(maxtime = 0),
#'    control_local = list(maxtime = 0))
#'    
#' hmm_biofam <- fit_biofam$model
#' }
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



