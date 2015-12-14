#' Three-channel biofam data
#'
#' Biofam data from the TraMineR package converted into three channels.
#'
#' @format A list including three sequence data sets for 2000 individuals with 16 state
#' variables, and a separate data frame with 1 id variable, 8 covariates, and 2 weights variables.
#'
#' @details This data is constructed from the \code{\link[TraMineR]{biofam}} data in
#' the TraMineR package. Here the original state sequences are converted into three separate
#' data sets: \code{children}, \code{married}, and \code{left}. These include the corresponding
#' life states from age 15 to 30: \code{childless} or (having) \code{children};
#' \code{single}, \code{married}, or \code{divorced}; and (living) \code{with parents} or
#' \code{left home}.
#'
#' Note that the \code{divorced} state does not give information on parenthood or residence,
#' so a guess is made based on preceeding states.
#'
#' The fourth data frame \code{covariates} is a collection of
#' additional variables from the original data:
#' \tabular{ll}{
#'  \code{idhous }\tab id\cr
#'  \code{sex }\tab sex\cr
#'  \code{birthyr }\tab birth year\cr
#'  \code{nat_1_02 }\tab first nationality\cr
#'  \code{plingu02 }\tab language of questionnaire\cr
#'  \code{p02r01 }\tab religion\cr
#'  \code{p02r04 }\tab religious participation\cr
#'  \code{cspfaj }\tab father's social status\cr
#'  \code{cspmoj }\tab mother's social status\cr
#'  \code{wp00tbgp }\tab weights inflating to the Swiss population\cr
#'  \code{wp00tbgs }\tab weights respecting sample size
#'  }
#'
#' The data is loaded by calling \code{data(biofam3c)}. It was built using following code:
#' \preformatted{
#' data("biofam" , package = "TraMineR")
#' biofam3c <- with(biofam, {
#'
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <- bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6 | bf == 7
#'
#' children[children == TRUE] <- "children"
#' children[children == FALSE] <- "childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 5) > 0) |
#'             (rowSums(bf == 7) > 0 & rowSums(bf == 6) > 0),]
#' children[rownames(bf) %in% rownames(div) & bf == 7] <- "children"
#'
#' married[married == TRUE] <- "married"
#' married[married == FALSE] <- "single"
#' married[bf == 7] <- "divorced"
#'
#' left[left == TRUE] <- "left home"
#' left[left == FALSE] <- "with parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 2) > 0 &
#'           rowSums(bf == 3) == 0 & rowSums(bf == 5) == 0 &
#'           rowSums(bf == 6) == 0) |
#'          (rowSums(bf == 7) > 0 & rowSums(bf == 4) > 0 &
#'           rowSums(bf == 3) == 0 & rowSums(bf == 5) == 0 &
#'           rowSums(bf == 6) == 0), ]
#' left[rownames(bf) %in% rownames(wp) & bf == 7] <- "with parents"
#'
#' list("children" = children, "married" = married, "left" = left,
#'   "covariates" = biofam[, c(1:9, 26:27)])
#' })
#' }
#'
#' @source \code{\link[TraMineR]{biofam}} data constructed from the Swiss Household Panel
#' \url{www.swisspanel.ch}
#'
#' @references Müller, N. S., M. Studer, G. Ritschard (2007). Classification de
#' parcours de vie à l'aide de l'optimal matching. In \emph{XIVe Rencontre de la Société
#' francophone de classification (SFC 2007), Paris, 5 - 7 septembre 2007}, pp. 157–160.
#'
#' @docType data
#' @keywords datasets
#' @name biofam3c
NULL
