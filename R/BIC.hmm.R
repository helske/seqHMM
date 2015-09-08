#' Compute Bayesian information criterion for a hidden Markov model or 
#' a mixture hidden Markov model
#' 
#' @export
#' @importFrom stats BIC
#' @rdname BIC
#' @param object Object of class \code{hmm} or \code{mhmm}.
#' @param trim Should trimming be performed before computing BIC? See
#'   \code{\link{trim_hmm}} for details. The default is \code{FALSE}.
#' @param ... further parameters to \code{\link{trim_hmm}}.
#'   
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and
#'   fitting Hidden Markov models, \code{\link{build_mhmm}} and 
#'   \code{\link{fit_mhmm}} for building and fitting mixture hidden Markov 
#'   models, and \code{\link{trim_hmm}} for finding better models by changing 
#'   small parameter values to zero.
#' 
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' # Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' 
#' # Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Starting values for emission matrices
#' B_child <- matrix(NA, nrow = 3, ncol = 2)
#' B_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 11:15])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow = 3, ncol = 2)
#' B_marr[1,] <- seqstatf(marr.seq[, 1:5])[, 2] + 0.1
#' B_marr[2,] <- seqstatf(marr.seq[, 6:10])[, 2] + 0.1
#' B_marr[3,] <- seqstatf(marr.seq[, 11:15])[, 2] + 0.1
#' B_marr <- B_marr / rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow = 3, ncol = 2)
#' B_left[1,] <- seqstatf(left.seq[, 1:5])[, 2] + 0.1
#' B_left[2,] <- seqstatf(left.seq[, 6:10])[, 2] + 0.1
#' B_left[3,] <- seqstatf(left.seq[, 11:15])[, 2] + 0.1
#' B_left <- B_left / rowSums(B_left)
#' 
#' # Starting values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Starting values for initial state probabilities
#' init <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = init
#'   )
#' 
#' # Fitting hidden Markov model
#' HMM <- fit_hmm(bHMM)
#' 
#' # Computing the BIC of the model
#' BIC(HMM$model)
#' 
BIC.hmm <- function(object, trim = FALSE,...){
  
  
  if(trim){
    trimmed <- trim_hmm(object, return_loglik=TRUE,...)
    object <- trimmed$model
    loglik <- trimmed$loglik
  } else loglik <- logLik(object)
  
  
  if(object$n_channels == 1){
    -2*loglik + log(object$n_sequences*object$length_of_sequences)*
      (sum(object$initial_probs>0)+sum(object$transition_matrix>0)+
         sum(object$emission_matrix>0))
  } else{  
    -2*loglik + log(object$n_sequences*object$length_of_sequences)*
      (sum(object$initial_probs>0)+sum(object$transition_matrix>0)+
         sum(sapply(object$emission_matrix,function(x) sum(x>0))))
  }
}