#' Three-channel biofam data
#'
#' Biofam data from the TraMineR package converted into three channels.
#'
#' @format A list including three sequence data sets for 2000 individuals with
#' 16 state variables, and a separate data frame with 1 id variable,
#' 8 covariates, and 2 weight variables.
#'
#' @details This data is constructed from the [TraMineR::biofam()]
#' data in the TraMineR package. Here the original state sequences are
#' converted into three separate data sets: `children`, `married`,
#' and `left`. These include the corresponding life states from age 15 to
#' 30: `childless` or (having) `children`; `single`,
#' `married`, or `divorced`; and (living) `with parents` or
#' `left home`.
#'
#' Note that the `divorced` state does not give information on parenthood
#' or residence, so a guess is made based on preceeding states.
#'
#' The fourth data frame `covariates` is a collection of
#' additional variables from the original data:
#' \tabular{ll}{
#'  `idhous `\tab id\cr
#'  `sex `\tab sex\cr
#'  `birthyr `\tab birth year\cr
#'  `nat_1_02 `\tab first nationality\cr
#'  `plingu02 `\tab language of questionnaire\cr
#'  `p02r01 `\tab religion\cr
#'  `p02r04 `\tab religious participation\cr
#'  `cspfaj `\tab father's social status\cr
#'  `cspmoj `\tab mother's social status\cr
#'  `wp00tbgp `\tab weights inflating to the Swiss population\cr
#'  `wp00tbgs `\tab weights respecting sample size
#'  }
#'
#' The data is loaded by calling `data(biofam3c)`. It was built using
#' following code:
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
#' children[rownames(bf) \%in\% rownames(div) & bf == 7] <- "children"
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
#' left[rownames(bf) \%in\% rownames(wp) & bf == 7] <- "with parents"
#'
#' list("children" = children, "married" = married, "left" = left,
#'   "covariates" = biofam[, c(1:9, 26:27)])
#' })
#' }
#'
#' @source [TraMineR::biofam()] data constructed from the Swiss
#' Household Panel
#' <https://forscenter.ch/projects/swiss-household-panel/>
#'
#' @references Müller, N. S., M. Studer, G. Ritschard (2007). Classification de
#' parcours de vie à l'aide de l'optimal matching. In *XIVe Rencontre de l
#' a Société francophone de classification (SFC 2007), Paris, 5 - 7 septembre
#' 2007*, pp. 157–160.
#'
#' @docType data
#' @keywords datasets
#' @name biofam3c
NULL
