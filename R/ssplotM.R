ssplotM <- function(x, hidden.paths=NULL, 
                   plots="obs", type="d", 
                   sortv=NULL, sort.channel=1, dist.method="OM",
                   with.missing=FALSE,
                   title=NA, title.n=TRUE, cex.title=1, title.pos=1,
                   withlegend="auto", ncol.legend="auto", 
                   with.missing.legend="auto",                         
                   legend.prop=0.3, cex.legend=1,
                   hidden.states.colors="auto", hidden.states.labels="auto",
                   xaxis=TRUE, xlab=NA, xtlab=NULL, xlab.pos=1,
                   ylab="auto", hidden.states.title="Hidden states", 
                   ylab.pos="auto", 
                   cex.lab=1, cex.axis=1, ...){
  # plot.new()
  grid.newpage()
  arguments <- as.list(match.call())[-1]
  if (length(arguments$x) == 1) {
    arguments$x <- arguments$x[[1]]
  }
  do.call(SSPlotter, args = do.call(ssp, args = arguments))
}