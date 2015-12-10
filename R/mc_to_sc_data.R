#' Merge Multiple Sequence Objects to One (Multichannel to Single Channel Data)
#'
#' Function \code{mc_to_sc_data} combines observed states of multiple sequence objects into one.
#'
#' @export
#' @param data A list of state sequence objects created with the \code{\link{seqdef}} function.
#' @param combine_missing Controls whether combined states of observations are
#'   coded missing (\code{NA}) if some of the channels include missing information.
#'   Defaults to \code{TRUE}.
#' @param all_combinations Controls whether all possible combinations of
#'   observed states are included in the single channel representation or only
#'   combinations that are found in the data. Defaults to \code{FALSE}, i.e.
#'   only actual observations are included.
#'
#' @examples
#'
#' ssplot(hmm_biofam$observations)
#' ssplot(mc_to_sc_data(hmm_biofam$observations))
#'
#' @seealso \code{\link{seqdef}} for creating state sequence objects.

mc_to_sc_data <- function(data, combine_missing = TRUE, all_combinations = FALSE){
  
  if (length(unique(sapply(data, nrow))) > 1) {
    stop("The number of subjects (rows) is not the same in all channels.")
  }
  if (length(unique(sapply(data, ncol))) > 1) {
    stop("The length of the sequences (number of columns) is not the same in all channels.")
  }

  alph <- apply(expand.grid(lapply(data,alphabet)), 1, paste0, collapse = "/")

  datax <- data[[1]]
  for(i in 2:length(data))
    datax <- as.data.frame(mapply(paste, datax,
      data[[i]], USE.NAMES=FALSE,SIMPLIFY=FALSE,
      MoreArgs=list(sep="/")))
  names(datax) <- names(data[[1]])
  if(combine_missing==TRUE){
    datax[Reduce("|", lapply(data,
        function(x)
          x==attr(data[[1]], "nr") |
          x==attr(data[[1]], "void") |
          is.na(x)))]<-NA
  }

  if (length(alph) <= 200) {
    cpal <- seqHMM::colorpalette[[length(alph)]]
  } else {
    cp <- NULL
    k <- 200
    p <- 0
    while(length(alph) - p > 0){
      cp <- c(cp, seqHMM::colorpalette[[k]])
      p <- p + k
      k <- k - 1
    }
    cpal <- cp[1:length(alph)]
  }


  if(all_combinations==TRUE){
    datax <- suppressWarnings(suppressMessages(seqdef(datax, alphabet=alph)))
  }else{
    datax<-suppressWarnings(suppressMessages((seqdef(datax))))
  }


  attr(datax, "xtstep") <- attr(data[[1]], "xtstep")
  attr(datax, "missing.color") <- attr(data[[1]], "missing.color")
  attr(datax, "nr") <- attr(data[[1]], "nr")
  attr(datax, "void") <- attr(data[[1]], "void")
  attr(datax, "missing") <- attr(data[[1]], "missing")
  attr(datax, "start") <- attr(data[[1]], "start")
  attr(datax, "cpal") <- cpal[alph %in% alphabet(datax)]

  datax
}
