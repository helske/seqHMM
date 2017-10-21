#' Merge Multiple Sequence Objects into One (from Multichannel to Single Channel Data)
#'
#' Function \code{mc_to_sc_data} combines observed states of multiple 
#'   sequence objects into one, time point by time point.
#'
#' @export
#' @param data A list of state sequence objects (\code{stslist}s) 
#'   created with the \code{\link{seqdef}} function.
#' @param combine_missing Controls whether combined states of observations 
#'   at time t are coded missing (coded with * in \code{stslist}s) 
#'   if one or more of the channels include missing information at time t. 
#'   Defaults to \code{TRUE}. \code{FALSE} keeps missing states
#'   as they are, producing more states in data; e.g. single/childless/* 
#'   where the observation in channel 3 is missing. 
#' @param all_combinations Controls whether all possible combinations of
#'   observed states are included in the single channel representation or 
#'   only combinations that are found in the data. Defaults to \code{FALSE}, 
#'   i.e. only actual observations are included.
#'
#' @seealso \code{\link{mc_to_sc}} for transforming multichannel \code{hmm} 
#'   or \code{mhmm} objects into single-channel representations; 
#'   \code{\link{ssplot}} for plotting multiple sequence data sets in the
#'   same plot; and \code{\link{seqdef}} for creating state sequence objects.
#' 
#' @examples
#' # Load three-channel sequence data
#' data("biofam3c")
#' 
#' # Building sequence objects
#' marr_seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child_seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left_seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' # Define colors
#' attr(marr_seq, "cpal") <- c("violetred2", "darkgoldenrod2", "darkmagenta")
#' attr(child_seq, "cpal") <- c("darkseagreen1", "coral3")
#' attr(left_seq, "cpal") <- c("lightblue", "red3")
#'
#' # Converting multichannel data to single-channel data
#' sc_data <- mc_to_sc_data(list(marr_seq, child_seq, left_seq))
#' 
#' # 10 combined states
#' alphabet(sc_data)
#' 
#' # Colors for combined states
#' attr(sc_data, "cpal") <- colorpalette[[14]][1:10]
#'
#' # Plotting sequences for the first 10 subjects
#' ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq, 
#'   "Residence" = left_seq, "Combined" = sc_data), type = "I",
#'   tlim = 1:10)
#' 
#' 
#' # Including all combinations (whether or not available in data)
#' sc_data_all <- mc_to_sc_data(list(marr_seq, child_seq, left_seq),
#'   all_combinations = TRUE)
#' 
#' # 12 combined states, 2 with no observations in data
#' seqstatf(sc_data_all)
#' 

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
