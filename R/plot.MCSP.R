#' Plot Multichannel Sequences and/or Most Probable Paths from Hidden Markov 
#' Models
#' 
#' Function \code{plot.MCSP} plots MCSP objects defined with
#' \code{\link{defineMCSP}}
#' 
#' 
#' 
#' @export
#' @param x A MCSP object.
#' @param ... Ignored.
#'   
#' @seealso \code{\link{defineMCSP}} for defining the plot before using 
#'   \code{plot.MCSP}, \code{\link{gridplot}} for plotting multiple MCSP objects.
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
#' ## Choosing colors
#' attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(marr.seq, "cpal") <- c("#E7298A", "#E6AB02")
#' attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#' 
#' 
#' # Plotting state distribution plots of observations
#' plot(defineMCSP(list(child.seq, marr.seq, left.seq), type="d", plots="obs"))


plot.MCSP <- function(x, ...) {
  plot.new()  
  grid.newpage()
  savepar <- par(no.readonly = TRUE)
  do.call(MCSPlotter,args=x)
  par(savepar)
}


