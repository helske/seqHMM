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
#' The data is created by calling \code{data(biofam3c)} using the following code:
#' 
#' \preformatted{
#' require(TraMineR)
#' data(biofam)
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf==4 | bf==5 | bf==6
#' married <- bf == 2 | bf== 3 | bf==6
#' left <- bf==1 | bf==3 | bf==5 | bf==6 | bf==7
#' 
#' children[children==TRUE] <- "children"
#' children[children==FALSE] <- "childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf==7)>0 & rowSums(bf==5)>0) |
#'             (rowSums(bf==7)>0 & rowSums(bf==6)>0),]
#' children[rownames(bf) %in% rownames(div) & bf==7] <- "children"
#' 
#' married[married==TRUE] <- "married"
#' married[married==FALSE] <- "single"
#' married[bf==7] <- "divorced"
#' 
#' left[left==TRUE] <- "left home"
#' left[left==FALSE] <- "with parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf==7)>0 & rowSums(bf==2)>0 & rowSums(bf==3)==0 &
#'             rowSums(bf==5)==0 & rowSums(bf==6)==0) |
#'            (rowSums(bf==7)>0 & rowSums(bf==4)>0 & rowSums(bf==3)==0 &
#'               rowSums(bf==5)==0 & rowSums(bf==6)==0),]
#' left[rownames(bf) %in% rownames(wp) & bf==7] <- "with parents"
#' 
#' biofam3c <- list("children" = children, "married" = married, 
#'                  "left" = left, "covariates" = biofam[, c(1:9, 26:27)])
#' }
#' 
#' @source \code{\link{biofam}} data constructed from the Swiss Household Panel 
#' \url{www.swisspanel.ch}
#' 
#' @references M{\"u}ller, N. S., M. Studer, G. Ritschard (2007). Classification de 
#' parcours de vie à l'aide de l'optimal matching. In \emph{XIVe Rencontre de la Société 
#' francophone de classification (SFC 2007), Paris, 5 - 7 septembre 2007}, pp. 157–160.
#' 
#' @docType data
#' @keywords datasets
#' @name biofam3c