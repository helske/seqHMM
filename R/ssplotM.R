ssplotM <- function(x, mpp=NULL, 
                   plots="obs", type="I", 
                   sortv=NULL, sort.channel=1, dist.method="OM",
                   with.missing=FALSE,
                   title=NA, title.n=TRUE, cex.title=1, title.pos=1,
                   withlegend="auto", ncol.legend="auto", 
                   with.missing.legend="auto",                         
                   legend.prop=0.3, cex.legend=1,
                   mpp.color="auto", mpp.labels="auto",
                   xaxis=TRUE, xlab=NA, xtlab=NULL, xlab.pos=1,
                   ylab="auto", hiddenStates.title="Hidden states", 
                   ylab.pos="auto", 
                   cex.lab=1, cex.axis=1, ...){
  # plot.new()
  grid.newpage()
  # sspargs <- do.call(ssp,args=as.list(match.call())[-1])
  do.call(SSPlotter,args=sspargs <- do.call(ssp,args=as.list(match.call())[-1]))
}