#' Mixture hidden Markov model for the biofam data
#' 
#' A mixture hidden Markov model (MHMM) fitted for the \code{\link[TraMineR]{biofam}} data.
#' 
#' @format A mixture hidden Markov model of class \code{mhmm}: 
#' three clusters with left-to-right models including 4, 4, and 6 hidden states.
#' Two covariates, \code{sex} and \code{cohort}, explaining cluster membership.
#'   
#' 
#' @details 
#' The model is loaded by calling \code{data(mhmm_biofam)}. It was created with the 
#' following code:
#' \preformatted{
#' require(TraMineR)
#' data(biofam3c)
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(biofam3c$children, start = 15)
#' marr.seq <- seqdef(biofam3c$married, start = 15)
#' left.seq <- seqdef(biofam3c$left, start = 15)
#' 
#' ## Choosing colors
#' attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#' 
#' ## Starting values for emission probabilities
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' emiss_1_child <- matrix(
#'   c(0.99, 0.01, # High probability for childless
#'     0.99, 0.01,
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' alphabet(marr.seq)
#' emiss_1_marr <- matrix(
#'   c(0.01, 0.01, 0.98, # High probability for single
#'     0.01, 0.01, 0.98,
#'     0.01, 0.98, 0.01, # High probability for married
#'     0.98, 0.01, 0.01), # High probability for divorced
#'   nrow = 4, ncol = 3, byrow = TRUE)
#' 
#' alphabet(left.seq)
#' emiss_1_left <- matrix(
#'   c(0.01, 0.99, # High probability for living with parents
#'     0.99, 0.01, # High probability for having left home
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' # Cluster 2
#' emiss_2_child <- matrix(
#'   c(0.99, 0.01, # High probability for childless
#'     0.99, 0.01,
#'     0.99, 0.01,
#'     0.01, 0.99), 
#'   nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' emiss_2_marr <- matrix(
#'   c(0.01, 0.01, 0.98, # High probability for single
#'     0.01, 0.01, 0.98,
#'     0.01, 0.98, 0.01, # High probability for married
#'     0.29, 0.7, 0.01),
#'   nrow = 4, ncol = 3, byrow = TRUE)
#' 
#' emiss_2_left <- matrix(
#'   c(0.01, 0.99, # High probability for living with parents
#'     0.99, 0.01,
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' # Cluster 3
#' emiss_3_child <- matrix(
#'   c(0.99, 0.01, # High probability for childless
#'     0.99, 0.01,
#'     0.01, 0.99,
#'     0.99, 0.01,
#'     0.01, 0.99,
#'     0.01, 0.99), 
#'   nrow = 6, ncol = 2, byrow = TRUE)
#' 
#' emiss_3_marr <- matrix(
#'   c(0.01, 0.01, 0.98, # High probability for single
#'     0.01, 0.01, 0.98,
#'     0.01, 0.01, 0.98,
#'     0.01, 0.98, 0.01,
#'     0.01, 0.98, 0.01, # High probability for married
#'     0.98, 0.01, 0.01), # High probability for divorced
#'   nrow = 6, ncol = 3, byrow = TRUE)
#' 
#' emiss_3_left <- matrix(
#'   c(0.01, 0.99, # High probability for living with parents
#'     0.99, 0.01,
#'     0.50, 0.50,
#'     0.01, 0.99,
#'     0.99, 0.01,
#'     0.99, 0.01), 
#'   nrow = 6, ncol = 2, byrow = TRUE)
#' 
#' # Initial values for transition matrices
#' trans_1 <- matrix(
#'   c(0.80, 0.16, 0.03, 0.01,
#'        0, 0.90, 0.07, 0.03,
#'        0,    0, 0.90, 0.10,
#'        0,    0,    0,    1),
#'   nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' trans_2 <- matrix(
#'   c(0.80, 0.10, 0.05,  0.03, 0.01, 0.01,
#'        0, 0.70, 0.10,  0.10, 0.05, 0.05,
#'        0,    0, 0.85,  0.01, 0.10, 0.04,
#'        0,    0,    0,  0.90, 0.05, 0.05,
#'        0,    0,    0,     0, 0.90,  0.1,
#'        0,    0,    0,     0,    0,    1),
#'   nrow = 6, ncol = 6, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities
#' initialProbs_1 <- c(0.9, 0.07, 0.02, 0.01)
#' initialProbs_2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Birth cohort
#' biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
#' biofam3c$covariates$cohort <- factor(
#'   biofam3c$covariates$cohort, labels=c("1909-1935", "1936-1945", "1946-1957"))
#' 
#' # Build mixture HMM
#' init_mhmm_biofam <- build_mhmm(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_probs = list(trans_1,trans_1,trans_2),
#'   emission_probs = list(list(emiss_1_child, emiss_1_marr, emiss_1_left),
#'                        list(emiss_2_child, emiss_2_marr, emiss_2_left),
#'                        list(emiss_3_child, emiss_3_marr, emiss_3_left)),
#'   initial_probs = list(initialProbs_1, initialProbs_1, initialProbs_2),
#'   formula = ~sex + cohort, data = biofam3c$covariates,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Parenthood", "Marriage", "Residence"),
#'   state_names = list(paste("State", 1:4), paste("State", 1:4), 
#'                      paste("State", 1:6)))
#' 
#' # Fitting the model
#' fit_biofam <- fit_mhmm(init_mhmm_biofam)
#' 
#' # Trimming the model
#' mhmm_biofam <- trim_hmm(fit_biofam$model, zerotol = 1e-04)
#' }
#'   
#' @seealso Examples of building and fitting MHMMs in \code{\link{build_mhmm}} and 
#' \code{\link{fit_mhmm}}; and \code{\link[TraMineR]{biofam}} for the original data and
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
#' \dontrun{
#' # Plotting the model for each cluster (change with Enter)
#' plot(mhmm_biofam)
#' }
#'   
NULL



