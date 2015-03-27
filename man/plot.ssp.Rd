% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/plot.ssp.R
\name{plot.ssp}
\alias{plot.ssp}
\title{Stack Multichannel Sequence Plots and/or Most Probable Paths Plots from Hidden Markov
Models}
\usage{
\method{plot}{ssp}(x, ...)
}
\arguments{
\item{x}{A ssp object.}

\item{...}{Ignored.}
}
\description{
Function \code{plot.ssp} plots stacked sequence plots from ssp objects defined with
\code{\link{ssp}}.
}
\examples{
require(TraMineR)

data(biofam)
biofam <- biofam[1:500,]

## Building one channel per type of event left, children or married
bf <- as.matrix(biofam[, 10:25])
children <-  bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6

children[children==TRUE] <- "Children"
children[children==FALSE] <- "Childless"

married[married==TRUE] <- "Married"
married[married==FALSE] <- "Single"

left[left==TRUE] <- "Left home"
left[left==FALSE] <- "With parents"

## Building sequence objects
child.seq <- seqdef(children)
marr.seq <- seqdef(married)
left.seq <- seqdef(left)

## Choosing colors
attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
attr(marr.seq, "cpal") <- c("#E7298A", "#E6AB02")
attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")


# Plotting state distribution plots of observations
ssp1 <- ssp(list(child.seq, marr.seq, left.seq), type="d", plots="obs")
plot(ssp1)
}
\seealso{
\code{\link{ssp}} for defining the plot before using
  \code{plot.ssp}, \code{\link{ssplot}} for straight plotting of sspobjects, \code{\link{gridplot}} for plotting multiple ssp objects.
}
