#' A feedback-augmented non-homogeneous hidden Markov Model for leaves data
#'
#' A FAN-HMM fitted for the`leaes` data.
#'
#' @format A model of class `fanhmm` with three hidden states
#'
#' @details
#' The model is loaded by calling `data(fanhmm_leaves)`. The code used to 
#' estimate the model is available on Github in `data-raw` folder.
#' @docType data
#' @keywords datasets
#' @name fanhmm_leaves
#' @examples
#' data("fanhmm_leaves")
#'
#' fanhmm_leaves
#' 
#' get_marginals(fanhmm_leaves)
#'
NULL
#' Synthetic data on fathers' parental leaves in Finland
#' 
#' @details The `leaves` data is a synthetic version of the Finnish fathers' leave-taking
#' data used in Helske et al. (2024) and Helske (2025). The data consists of 
#' variables
#'  * workplace: Workplace ID.
#'  * father: father ID within workplace. More accurately, this is the birth of 
#'  a child, i.e. same father can have multiple entries in data, but each entry 
#'  has separate ID.
#'  * year: Year when the child was born.
#'  * leave: Factor of leave-taking of the father.
#'  * Occupation: Factor of skill level of the father's occupation
#'  * reform2013: Factor indicating whether the father was eligible for the leave 
#'  under the 2013 reform.
#'  * same_occupation: Logical value, TRUE if father had same occupation as the 
#'  previous father.
#'  * lag_reform2013: Factor indicating whether the previous father was 
#'  eligible for the reform.
#'  * lag_occupation: Factor indiciting the occupation of previous father. 
#' 
#' @format A `data.table` with 9281 rows and 9 variables
#' @references
#' Helske S, Helske J, Chapman SN, Kotimäki S, Salin M, and Tikka S (2024). 
#' Heterogeneous workplace peer effects in fathers’ parental leave uptake in Finland. 
#' doi: 10.31235/osf.io/p3chf
#' Helske J (2025). Feedback-augmented Non-homogeneous Hidden Markov Models for Longitudinal Causal Inference.
#' ArXiv preprint. doi:10.48550/arXiv.2503.16014
#' @docType data
#' @keywords datasets
#' @name leaves
#' @examples
#' data("leaves")
#' head(leaves)
#' # convert to stslist
#' leaves_sequences <- data_to_stslist(
#'   leaves, id = "workplace", time = "father", responses = "leave",
#'   seqdef_args = list(cpal = c("tomato", "navyblue", "goldenrod"))
#' )
#' stacked_sequence_plot(leaves_sequences)
NULL
