ssplotM <- function(x, hidden.paths, 
                   plots, type, tlim,
                   sortv, sort.channel, dist.method,
                   with.missing,
                   title, title.n, cex.title, title.pos,
                   withlegend, ncol.legend, 
                   with.missing.legend,                         
                   legend.prop, cex.legend,
                   hidden.states.colors, hidden.states.labels,
                   xaxis, xlab, xtlab, xlab.pos,
                   yaxis, ylab, hidden.states.title, 
                   ylab.pos, 
                   cex.lab, cex.axis, ...){
  # plot.new()
  grid.newpage()
  arguments <- as.list(match.call())[-1]
  if (length(arguments$x) == 1) {
    arguments$x <- arguments$x[[1]]
  }
  do.call(SSPlotter, args = do.call(ssp, args = arguments))
}