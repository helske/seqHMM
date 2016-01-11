#' Mixture hidden Markov model for the biofam data
#'
#' A mixture hidden Markov model (MHMM) fitted for the \code{\link[TraMineR]{biofam}} data.
#'
#' @format A mixture hidden Markov model of class \code{mhmm}:
#' three clusters with left-to-right models including 4, 4, and 6 hidden states.
#' Two covariates, \code{sex} and \code{cohort}, explaining the cluster membership.
#'
#'
#' @details
#' The model was created with the  following code:
#' \preformatted{
#' data("biofam3c")
#'
#' ## Building sequence objects
#' marr_seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child_seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left_seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' ## Choosing colors
#' attr(marr_seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(child_seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(left_seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#'
#' ## Starting values for emission probabilities
#' # Cluster 1
#' B1_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.3, 0.6, 0.1, # High probability for married
#'     0.3, 0.3, 0.4), # High probability for divorced
#'   nrow = 4, ncol = 3, byrow = TRUE)
#'
#' B1_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.9, 0.1),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' B1_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9, # High probability for having left home
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' # Cluster 2
#'
#' B2_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.7, 0.2, 0.1),
#'   nrow = 4, ncol = 3, byrow = TRUE)
#'
#' B2_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' B2_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' # Cluster 3
#' B3_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.3, 0.4, 0.3,
#'     0.1, 0.1, 0.8), # High probability for divorced
#'   nrow = 6, ncol = 3, byrow = TRUE)
#'
#' B3_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9),
#'   nrow = 6, ncol = 2, byrow = TRUE)
#'
#'
#' B3_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 6, ncol = 2, byrow = TRUE)
#'
#' # Starting values for transition matrices
#' A1 <- matrix(
#'   c(0.80, 0.16, 0.03, 0.01,
#'     0,    0.90, 0.07, 0.03,
#'     0,    0,    0.90, 0.10,
#'     0,    0,    0,       1),
#'   nrow = 4, ncol = 4, byrow = TRUE)
#'
#' A2 <- matrix(
#'   c(0.80, 0.10, 0.05, 0.03, 0.01, 0.01,
#'     0,    0.70, 0.10, 0.10, 0.05, 0.05,
#'     0,    0,    0.85, 0.01, 0.10, 0.04,
#'     0,    0,    0,    0.90, 0.05, 0.05,
#'     0,    0,    0,    0,    0.90, 0.10,
#'     0,    0,    0,    0,    0,       1),
#'   nrow = 6, ncol = 6, byrow = TRUE)
#'
#' # Starting values for initial state probabilities
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#'
#' # Birth cohort
#' biofam3c$covariates$cohort <- factor(cut(biofam3c$covariates$birthyr,
#'     c(1908, 1935, 1945, 1957)), labels = c("1909-1935", "1936-1945", "1946-1957"))
#'
#' # Build mixture HMM
#' init_mhmm_bf <- build_mhmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   transition_probs = list(A1, A1, A2),
#'   emission_probs = list(list(B1_marr, B1_child, B1_left),
#'     list(B2_marr, B2_child, B2_left),
#'     list(B3_marr, B3_child, B3_left)),
#'   formula = ~sex + cohort, data = biofam3c$covariates,
#'   channel_names = c("Marriage", "Parenthood", "Residence"))
#'
#' # Fitting the model
#' mhmm_biofam <- fit_model(init_mhmm_bf)$model
#' }
#'
#' @seealso Examples of building and fitting MHMMs in \code{\link{build_mhmm}} and
#' \code{\link{fit_model}}; and \code{\link[TraMineR]{biofam}} for the original data and
#' \code{\link{biofam3c}} for the three-channel version used in this model.
#'
#' @docType data
#' @keywords datasets
#' @name mhmm_biofam
#' @examples
#' data("mhmm_biofam")
#'
#' # use conditional_se = FALSE for more accurate standard errors
#' # (these are considerebly slower to compute)
#' summary(mhmm_biofam$model)
#'
#' if (interactive()) {
#'   # Plotting the model for each cluster (change with Enter)
#'   plot(mhmm_biofam)
#' }
#'
NULL



