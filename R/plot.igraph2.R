# A copy of plot.igraph from igraph's development version
# We need to use internal functions of igraph as the version on CRAN does not yet
# contain a requested feature (available on GitHub)
plot.igraph2 <- function(x, axes = FALSE, add = FALSE, xlim = c(-1, 1),
                         ylim = c(-1, 1), mark.groups = list(), mark.shape = 1/2,
                         mark.col = rainbow(length(mark.groups), alpha = 0.3),
                         mark.border = rainbow(length(mark.groups), alpha = 1),
                         mark.expand = 15, ...) {
  graph <- x
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  params <- igraph:::i.parse.plot.params(graph, list(...))
  vertex.size <- 1/200 * params("vertex", "size")
  label.family <- params("vertex", "label.family")
  label.font <- params("vertex", "label.font")
  label.cex <- params("vertex", "label.cex")
  label.degree <- params("vertex", "label.degree")
  label.color <- params("vertex", "label.color")
  label.dist <- params("vertex", "label.dist")
  labels <- params("vertex", "label")
  shape <- igraph:::igraph.check.shapes(params("vertex", "shape"))
  edge.color <- params("edge", "color")
  edge.width <- params("edge", "width")
  edge.lty <- params("edge", "lty")
  arrow.mode <- params("edge", "arrow.mode")
  edge.labels <- params("edge", "label")
  loop.angle <- params("edge", "loop.angle")
  edge.label.font <- params("edge", "label.font")
  edge.label.family <- params("edge", "label.family")
  edge.label.cex <- params("edge", "label.cex")
  edge.label.color <- params("edge", "label.color")
  elab.x <- params("edge", "label.x")
  elab.y <- params("edge", "label.y")
  arrow.size <- params("edge", "arrow.size")[1]
  arrow.width <- params("edge", "arrow.width")[1]
  curved <- params("edge", "curved")
  if (is.function(curved)) {
    curved <- curved(graph)
  }
  layout <- params("plot", "layout")
  margin <- params("plot", "margin")
  margin <- rep(margin, length = 4)
  rescale <- params("plot", "rescale")
  asp <- params("plot", "asp")
  frame <- params("plot", "frame")
  main <- params("plot", "main")
  sub <- params("plot", "sub")
  xlab <- params("plot", "xlab")
  ylab <- params("plot", "ylab")
  palette <- params("plot", "palette")
  if (!is.null(palette)) {
    old_palette <- palette(palette)
    on.exit(palette(old_palette), add = TRUE)
  }
  arrow.mode <- igraph:::i.get.arrow.mode(graph, arrow.mode)
  maxv <- max(vertex.size)
  if (rescale) {
    layout <- norm_coords(layout, -1, 1, -1, 1)
    xlim <- c(xlim[1] - margin[2] - maxv, xlim[2] + margin[4] +
                maxv)
    ylim <- c(ylim[1] - margin[1] - maxv, ylim[2] + margin[3] +
                maxv)
  }
  if (!add) {
    plot(0, 0, type = "n", xlab = xlab, ylab = ylab, xlim = xlim,
         ylim = ylim, axes = axes, frame = frame, asp = asp,
         main = main, sub = sub)
  }
  if (!is.list(mark.groups) && is.numeric(mark.groups)) {
    mark.groups <- list(mark.groups)
  }
  mark.shape <- rep(mark.shape, length = length(mark.groups))
  mark.border <- rep(mark.border, length = length(mark.groups))
  mark.col <- rep(mark.col, length = length(mark.groups))
  mark.expand <- rep(mark.expand, length = length(mark.groups))
  for (g in seq_along(mark.groups)) {
    v <- V(graph)[mark.groups[[g]]]
    if (length(vertex.size) == 1) {
      vs <- vertex.size
    }
    else {
      vs <- rep(vertex.size, length = vcount(graph))[v]
    }
    igraph:::igraph.polygon(layout[v, , drop = FALSE], vertex.size = vs,
                   expand.by = mark.expand[g]/200, shape = mark.shape[g],
                   col = mark.col[g], border = mark.border[g])
  }
  el <- as_edgelist(graph, names = FALSE)
  loops.e <- which(el[, 1] == el[, 2])
  nonloops.e <- which(el[, 1] != el[, 2])
  loops.v <- el[, 1][loops.e]
  loop.labels <- edge.labels[loops.e]
  loop.labx <- if (is.null(elab.x)) {
    rep(NA, length(loops.e))
  }
  else {
    elab.x[loops.e]
  }
  loop.laby <- if (is.null(elab.y)) {
    rep(NA, length(loops.e))
  }
  else {
    elab.y[loops.e]
  }
  edge.labels <- edge.labels[nonloops.e]
  elab.x <- if (is.null(elab.x))
    NULL
  else elab.x[nonloops.e]
  elab.y <- if (is.null(elab.y))
    NULL
  else elab.y[nonloops.e]
  el <- el[nonloops.e, , drop = FALSE]
  edge.coords <- matrix(0, nrow = nrow(el), ncol = 4)
  edge.coords[, 1] <- layout[, 1][el[, 1]]
  edge.coords[, 2] <- layout[, 2][el[, 1]]
  edge.coords[, 3] <- layout[, 1][el[, 2]]
  edge.coords[, 4] <- layout[, 2][el[, 2]]
  if (length(unique(shape)) == 1) {
    ec <- igraph:::.igraph.shapes[[shape[1]]]$clip(edge.coords, el,
                                          params = params, end = "both")
  }
  else {
    shape <- rep(shape, length = vcount(graph))
    ec <- edge.coords
    ec[, 1:2] <- t(sapply(seq(length = nrow(el)), function(x) {
      igraph:::.igraph.shapes[[shape[el[x, 1]]]]$clip(edge.coords[x,
                                                         , drop = FALSE], el[x, , drop = FALSE], params = params,
                                             end = "from")
    }))
    ec[, 3:4] <- t(sapply(seq(length = nrow(el)), function(x) {
      igraph:::.igraph.shapes[[shape[el[x, 2]]]]$clip(edge.coords[x,
                                                         , drop = FALSE], el[x, , drop = FALSE], params = params,
                                             end = "to")
    }))
  }
  x0 <- ec[, 1]
  y0 <- ec[, 2]
  x1 <- ec[, 3]
  y1 <- ec[, 4]
  if (length(loops.e) > 0) {
    ec <- edge.color
    if (length(ec) > 1) {
      ec <- ec[loops.e]
    }
    point.on.cubic.bezier <- function(cp, t) {
      c <- 3 * (cp[2, ] - cp[1, ])
      b <- 3 * (cp[3, ] - cp[2, ]) - c
      a <- cp[4, ] - cp[1, ] - c - b
      t2 <- t * t
      t3 <- t * t * t
      a * t3 + b * t2 + c * t + cp[1, ]
    }
    compute.bezier <- function(cp, points) {
      dt <- seq(0, 1, by = 1/(points - 1))
      sapply(dt, function(t) point.on.cubic.bezier(cp,
                                                   t))
    }
    plot.bezier <- function(cp, points, color, width, arr,
                            lty, arrow.size, arr.w) {
      p <- compute.bezier(cp, points)
      polygon(p[1, ], p[2, ], border = color, lwd = width,
              lty = lty)
      if (arr == 1 || arr == 3) {
        igraph:::igraph.Arrows(p[1, ncol(p) - 1], p[2, ncol(p) -
                                             1], p[1, ncol(p)], p[2, ncol(p)], sh.col = color,
                      h.col = color, size = arrow.size, sh.lwd = width,
                      h.lwd = width, open = FALSE, code = 2, width = arr.w)
      }
      if (arr == 2 || arr == 3) {
        igraph:::igraph.Arrows(p[1, 2], p[2, 2], p[1, 1], p[2,
                                                   1], sh.col = color, h.col = color, size = arrow.size,
                      sh.lwd = width, h.lwd = width, open = FALSE,
                      code = 2, width = arr.w)
      }
    }
    loop <- function(x0, y0, cx = x0, cy = y0, color, angle = 0,
                     label = NA, width = 1, arr = 2, lty = 1, arrow.size = arrow.size,
                     arr.w = arr.w, lab.x, lab.y) {
      rad <- angle
      center <- c(cx, cy)
      cp <- matrix(c(x0, y0, x0 + 0.4, y0 + 0.2, x0 + 0.4,
                     y0 - 0.2, x0, y0), ncol = 2, byrow = TRUE)
      phi <- atan2(cp[, 2] - center[2], cp[, 1] - center[1])
      r <- sqrt((cp[, 1] - center[1])^2 + (cp[, 2] - center[2])^2)
      phi <- phi + rad
      cp[, 1] <- cx + r * cos(phi)
      cp[, 2] <- cy + r * sin(phi)
      plot.bezier(cp, 50, color, width, arr = arr, lty = lty,
                  arrow.size = arrow.size, arr.w = arr.w)
      if (is.language(label) || !is.na(label)) {
        lx <- x0 + 0.3
        ly <- y0
        phi <- atan2(ly - center[2], lx - center[1])
        r <- sqrt((lx - center[1])^2 + (ly - center[2])^2)
        phi <- phi + rad
        lx <- cx + r * cos(phi)
        ly <- cy + r * sin(phi)
        if (!is.na(lab.x)) {
          lx <- lab.x
        }
        if (!is.na(lab.y)) {
          ly <- lab.y
        }
        text(lx, ly, label, col = edge.label.color, font = edge.label.font,
             family = edge.label.family, cex = edge.label.cex)
      }
    }
    ec <- edge.color
    if (length(ec) > 1) {
      ec <- ec[loops.e]
    }
    vs <- vertex.size
    if (length(vertex.size) > 1) {
      vs <- vs[loops.v]
    }
    ew <- edge.width
    if (length(edge.width) > 1) {
      ew <- ew[loops.e]
    }
    la <- loop.angle
    if (length(loop.angle) > 1) {
      la <- la[loops.e]
    }
    lty <- edge.lty
    if (length(edge.lty) > 1) {
      lty <- lty[loops.e]
    }
    arr <- arrow.mode
    if (length(arrow.mode) > 1) {
      arr <- arrow.mode[loops.e]
    }
    asize <- arrow.size
    if (length(arrow.size) > 1) {
      asize <- arrow.size[loops.e]
    }
    xx0 <- layout[loops.v, 1] + cos(la) * vs
    yy0 <- layout[loops.v, 2] - sin(la) * vs
    mapply(loop, xx0, yy0, color = ec, angle = -la, label = loop.labels,
           lty = lty, width = ew, arr = arr, arrow.size = asize,
           arr.w = arrow.width, lab.x = loop.labx, lab.y = loop.laby)
  }
  if (length(x0) != 0) {
    if (length(edge.color) > 1) {
      edge.color <- edge.color[nonloops.e]
    }
    if (length(edge.width) > 1) {
      edge.width <- edge.width[nonloops.e]
    }
    if (length(edge.lty) > 1) {
      edge.lty <- edge.lty[nonloops.e]
    }
    if (length(arrow.mode) > 1) {
      arrow.mode <- arrow.mode[nonloops.e]
    }
    if (length(arrow.size) > 1) {
      arrow.size <- arrow.size[nonloops.e]
    }
    if (length(curved) > 1) {
      curved <- curved[nonloops.e]
    }
    if (length(unique(arrow.mode)) == 1) {
      lc <- igraph:::igraph.Arrows(x0, y0, x1, y1, h.col = edge.color,
                          sh.col = edge.color, sh.lwd = edge.width, h.lwd = 1,
                          open = FALSE, code = arrow.mode[1], sh.lty = edge.lty,
                          h.lty = 1, size = arrow.size, width = arrow.width,
                          curved = curved)
      lc.x <- lc$lab.x
      lc.y <- lc$lab.y
    }
    else {
      curved <- rep(curved, length = ecount(graph))[nonloops.e]
      lc.x <- lc.y <- numeric(length(curved))
      for (code in 0:3) {
        valid <- arrow.mode == code
        if (!any(valid)) {
          next
        }
        ec <- edge.color
        if (length(ec) > 1) {
          ec <- ec[valid]
        }
        ew <- edge.width
        if (length(ew) > 1) {
          ew <- ew[valid]
        }
        el <- edge.lty
        if (length(el) > 1) {
          el <- el[valid]
        }
        lc <- igraph:::igraph.Arrows(x0[valid], y0[valid], x1[valid],
                            y1[valid], code = code, sh.col = ec, h.col = ec,
                            sh.lwd = ew, h.lwd = 1, h.lty = 1, sh.lty = el,
                            open = FALSE, size = arrow.size, width = arrow.width,
                            curved = curved[valid])
        lc.x[valid] <- lc$lab.x
        lc.y[valid] <- lc$lab.y
      }
    }
    if (!is.null(elab.x)) {
      lc.x <- ifelse(is.na(elab.x), lc.x, elab.x)
    }
    if (!is.null(elab.y)) {
      lc.y <- ifelse(is.na(elab.y), lc.y, elab.y)
    }
    text(lc.x, lc.y, labels = edge.labels, col = edge.label.color,
         family = edge.label.family, font = edge.label.font,
         cex = edge.label.cex)
  }
  rm(x0, y0, x1, y1)
  if (length(unique(shape)) == 1) {
    igraph:::.igraph.shapes[[shape[1]]]$plot(layout, params = params)
  }
  else {
    sapply(seq(length = vcount(graph)), function(x) {
      igraph:::.igraph.shapes[[shape[x]]]$plot(layout[x, , drop = FALSE],
                                      v = x, params = params)
    })
  }
  par(xpd = TRUE)
  x <- layout[, 1] + label.dist * cos(-label.degree) * (vertex.size +
                                                          6 * 8)/200
  y <- layout[, 2] + label.dist * sin(-label.degree) * (vertex.size +
                                                          6 * 8)/200
  if (length(label.family) == 1) {
    text(x, y, labels = labels, col = label.color, family = label.family,
         font = label.font, cex = label.cex)
  }
  else {
    if1 <- function(vect, idx) if (length(vect) == 1)
      vect
    else vect[idx]
    sapply(seq_len(vcount(graph)), function(v) {
      text(x[v], y[v], labels = if1(labels, v), col = if1(label.color,
                                                          v), family = if1(label.family, v), font = if1(label.font,
                                                                                                        v), cex = if1(label.cex, v))
    })
  }
  rm(x, y)
  invisible(NULL)
}
