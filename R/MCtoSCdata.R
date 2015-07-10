#' Merge Multiple Sequence Objects to One (Multichannel to Single Channel Data)
#' 
#' Function \code{MCtoSCdata} combines observed states of multiple sequence objects into one.
#' 
#' @export
#' @param data A list of state sequence objects created with the \code{\link{seqdef}} function.
#' @param combine.missing Controls whether combined states of observations are 
#'   coded missing (NA) if some of the channels include missing information. 
#'   Defaults to \code{TRUE}.
#' @param all.combinations Controls whether all possible combinations of
#'   observed states are included in the single channel representation or only
#'   combinations that are found in the data. Defaults to \code{FALSE}, i.e.
#'   only actual observations are included.
#'   
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf==4 | bf==5 | bf==6
#' married <- bf == 2 | bf== 3 | bf==6
#' left <- bf==1 | bf==3 | bf==5 | bf==6
#' 
#' children[children==TRUE] <- "Children"
#' children[children==FALSE] <- "Childless"
#' 
#' married[married==TRUE] <- "Married"
#' married[married==FALSE] <- "Single"
#' 
#' left[left==TRUE] <- "Left home"
#' left[left==FALSE] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' scdata <- MCtoSCdata(list(child.seq, marr.seq, left.seq))
#' 
#' @seealso \code{\link{seqdef}} for creating state sequence objects.

MCtoSCdata <- function(data, combine.missing=TRUE, all.combinations=FALSE){

  alph <- apply(
    expand.grid(lapply(data,alphabet)),                
    1,paste0,collapse="/")
    
  datax <- data[[1]]
  for(i in 2:length(data))
    datax <- as.data.frame(mapply(paste, datax,
                                              data[[i]],
                                              USE.NAMES=FALSE,SIMPLIFY=FALSE,
                                              MoreArgs=list(sep="/")))
  names(datax) <- names(data[[1]])   
  if(combine.missing==TRUE){
    datax[Reduce("|",
                               lapply(
                                 data, 
                                 function(x) 
                                   x==attr(data[[1]], "nr") |
                                   x==attr(data[[1]], "void") |
                                   is.na(x)))]<-NA
  }
  
  cpal <- seqHMM::colorpalette[[length(alph)]]
  
  if(all.combinations==TRUE){
    datax <- suppressWarnings(suppressMessages(seqdef(datax, alphabet=alph)))
  }else{
    datax<-suppressWarnings(suppressMessages((seqdef(datax))))
  }
}