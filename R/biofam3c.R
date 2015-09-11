#' Three-channel biofam data
#' 
#' Biofam data from the TraMineR package converted into three channels.
#' 
#' @format A list including three sequence data sets for 2000 individuals with 16 state 
#' variables, and a separate data frame with 1 id variable, 8 covariates, and 2 weights variables.
#' 
#' @details This data is constructed from the \code{\link{biofam}} data in 
#' \code{\link{TraMineR}}. Here the original state sequences are converted into three separate
#' data sets: \code{children}, \code{married}, and \code{left}. These include the corresponding 
#' life states from age 15 to 30: \code{childless} or (having) \code{children}; 
#' \code{single}, \code{married}, or \code{divorced}; and (living) \code{with parents} or
#' \code{left home}. 
#' 
#' The fourth data frame \code{covariates} is a collection of 
#' additional variables from the original data:
#' 
#' \describe{
#'  \item{\code{idhous}}{id}
#'  \item{\code{sex}}{sex}
#'  \item{\code{birthyr}}{birth year}
#'  \item{\code{nat_1_02}}{first nationality}
#'  \item{\code{plingu02}}{language of questionnaire}
#'  \item{\code{p02r01}}{religion}
#'  \item{\code{p02r04}}{religious participation}
#'  \item{\code{cspfaj}}{father's social status}
#'  \item{\code{cspmoj}}{mother's social status}
#'  \item{\code{wp00tbgp}}{weights inflating to the Swiss population}
#'  \item{\code{wp00tbgs}}{weights respecting sample size}
#' }
#'   
#' @source \code{\link{biofam}} data constructed from the Swiss Household Panel 
#' \url{www.swisspanel.ch}
#'   
#' @seealso \code{\link{biofam}} for the original single-channel data available in 
#' the TraMineR package, \code{\link{seqdef}} for defining state sequence objects, and
#' \code{\link{ssplot}} for plotting multichannel sequence data.
#' @docType data
#' @keywords datasets
#' @name biofam3c