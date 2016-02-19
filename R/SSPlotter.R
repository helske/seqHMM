SSPlotter <- function(obs, nchannels, nplots,
  legend.c.prop, legend.r.prop, ylab.space, xaxis.space, xt.space,
  hidden.paths.seq, orderv, plotxaxis,
  hidden.paths, plots, type, tlim, n.seq, sortv, sort.channel,
  with.missing, title, title.n, cex.title, title.pos,
  withlegend, ncol.legend, with.missing.legend,
  legend.prop, cex.legend, hidden.states.colors, hidden.states.labels,
  xaxis, xlab, xtlab, xlab.pos, yaxis, ylab,
  hidden.states.title, ylab.pos,
  cex.lab, cex.axis, new, call = match.call(), ...){

  if (xlab.pos * cex.lab + xaxis.space + xt.space > 0) {
    xspace <- xlab.pos * cex.lab + xaxis.space + xt.space + 0.5
  } else {
    xspace <- 0
  }

  # Grid for plotting regions
  top.vp <- viewport(layout = grid.layout(nrow = 4, ncol = 3,
    widths = unit(c(ylab.space, 1, legend.c.prop),
      c("lines", "null", "npc")),
    heights = unit(c((title.pos*cex.title+1.5), 1,
      xspace,
      legend.r.prop),
      c("lines", "null", "lines", "npc"))))

  vptitle <- viewport(layout.pos.row = 1, layout.pos.col = 2, name = "vptitle")
  vpylab <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = "vpylab")
  vpplot <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = "vpplot")
  if (withlegend == "right" || withlegend == "right.combined") {
    vplegend <- viewport(layout.pos.row = 2, layout.pos.col = 3, name = "vplegend")
  } else if (withlegend == "bottom" || withlegend == "bottom.combined") {
    vplegend <- viewport(layout.pos.row = 4, layout.pos.col = 2, name = "vplegend")
  }

  vpxaxis <- viewport(layout.pos.row = 3, layout.pos.col = 2, name = "vpxaxis")

  if (withlegend == FALSE) {
    splot <- vpTree(top.vp, vpList(vptitle, vpylab, vpplot, vpxaxis))
  } else {
    splot <- vpTree(top.vp, vpList(vptitle, vpylab, vpplot, vplegend, vpxaxis))
  }

  pushViewport(splot)

  # Grid for plots
  upViewport()
  downViewport("vpplot")
  pushViewport(viewport(layout = grid.layout(nrow = nplots, ncol = 1),
    width = unit(0.95, "npc")))

  # Plotting observations
  if (plots == "both" || plots == "obs") {
    if (type == "I") {
      if (is.null(sortv)) {
        if (nchannels > 1) {
          for(i in 1:(nchannels-1)) {
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
            par(plt = gridPLT(), new = TRUE)
            seqplot(obs[[i]], type = "I", tlim = tlim, withlegend = FALSE,
              use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
              xtlab = xtlab, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[nchannels]], type = "I", tlim = tlim, withlegend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.plot = cex.axis, ...)
          popViewport()
        } else {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs, type = "I", tlim = tlim, withlegend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.plot = cex.axis, ...)
          popViewport()
        }

      } else if (length(sortv) == 1 && (sortv == "from.start" || sortv == "from.end")) {
        if (nchannels > 1) {
          for (i in 1:(nchannels-1)) {
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
            par(plt = gridPLT(), new = TRUE)
            seqplot(obs[[i]][orderv,], type = "I", tlim = tlim, withlegend = FALSE,
              use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
              xtlab = xtlab, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[nchannels]][orderv,], type = "I", tlim = tlim, withlegend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.plot = cex.axis, ...)
          popViewport()
        } else {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[orderv,], type = "I", tlim = tlim, withlegend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.plot = cex.axis, ...)
          popViewport()
        }


      } else if(length(sortv) > 1) {
        if (nchannels > 1) {
          for (i in 1:(nchannels-1)) {
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
            par(plt = gridPLT(), new = TRUE)
            seqplot(obs[[i]], type = "I", tlim = tlim, sortv = sortv, withlegend = FALSE,
              use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
              xtlab = xtlab, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[nchannels]], type = "I", tlim = tlim, sortv = sortv, withlegend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.plot = cex.axis, ...)
          popViewport()
        } else {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs, type = "I", tlim = tlim, sortv = sortv, withlegend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.plot = cex.axis, ...)
          popViewport()
        }

      }
      if(plots == "obs"){
        # Close viewport "vpplot"
        popViewport(2)
      }
    } else{  # if type == "d"
      if (nchannels > 1) {
        for (i in 1:(nchannels-1)) {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[i]], type = "d", withlegend = FALSE,
            use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
            with.missing = with.missing, xtlab = xtlab, ...)
          popViewport()
        }
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
        par(plt = gridPLT(), new = TRUE)
        seqplot(obs[[nchannels]], type = "d", withlegend = FALSE, xaxis = plotxaxis,
          use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
          with.missing = with.missing, xtlab = xtlab, cex.plot = cex.axis, ...)
        popViewport()
      } else {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
        par(plt = gridPLT(), new = TRUE)
        seqplot(obs, type = "d", withlegend = FALSE, xaxis = plotxaxis,
          use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
          with.missing = with.missing, xtlab = xtlab, cex.plot = cex.axis, ...)
        popViewport()
      }

      if(plots == "obs"){
        # Close viewport "vpplot"
        popViewport(2)
      }
    }
  }

  # Plotting the most probable paths
  if (plots == "both" || plots == "hidden.paths") {
    if (type == "I") {
      if (is.null(sortv)) {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
        par(plt = gridPLT(), new = TRUE)
        seqplot(hidden.paths.seq, type = type, tlim = tlim, sortv = sortv, withlegend = FALSE,
          use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
          xtlab = xtlab, cex.plot = cex.axis, ...)
        popViewport()
      } else if (length(sortv) == 1 && (sortv == "from.start" || sortv == "from.end")) {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
        par(plt = gridPLT(), new = TRUE)
        seqplot(hidden.paths.seq[orderv,], type = type, tlim = tlim, withlegend = FALSE,
          use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
          xtlab = xtlab, cex.plot = cex.axis, ...)
        popViewport()
      } else if(length(sortv) > 1) {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
        par(plt = gridPLT(), new = TRUE)
        seqplot(hidden.paths.seq, type = type, tlim = tlim, sortv = sortv, withlegend = FALSE,
          use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
          xtlab = xtlab, cex.plot = cex.axis, ...)
        popViewport()
      }
    } else {
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
      par(plt = gridPLT(), new = TRUE)
      seqplot(hidden.paths.seq, type = type, withlegend = FALSE,
        use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
        xtlab = xtlab, cex.plot = cex.axis, ...)
      popViewport()
    }
    # Close viewport "vpplot"
    popViewport(2)
  }

  # Plotting x label
  if (!is.na(xlab)) {
    downViewport("vpxaxis")
    grid.text(xlab, y = unit(1, "lines"),
      gp = gpar(cex = cex.lab))
    popViewport()
  }

  # Plotting y labels (channels and/or hidden states)
  if (length(ylab) > 1 || (length(ylab) == 1 && !is.na(ylab) && ylab != FALSE)) {
    downViewport("vpylab")
    pushViewport(viewport(layout = grid.layout(nrow = nplots, ncol = 1)))
    if(plots == "both" || plots == "obs"){
      for(i in 1:nchannels){
        pushViewport(viewport(layout.pos.row = i,
          layout.pos.col = 1))
        grid.text(ylab[i], x = unit(ylab.space/cex.lab-
            ylab.pos[i]+1, "lines"),
          gp = gpar(cex = cex.lab), rot = 90, vjust = 0.5, hjust = 0.5)
        popViewport()
      }
    }
    if (plots == "both" || plots == "hidden.paths") {
      pushViewport(viewport(layout.pos.row = nplots,
        layout.pos.col = 1))
      grid.text(hidden.states.title, x = unit(ylab.space/cex.lab-
          ylab.pos[nplots]+1,
        "lines"),
        gp = gpar(cex = cex.lab), rot = 90, vjust = 0.5, hjust = 0.5)
      popViewport()
    }
    # Close viewport "vpylab"
    popViewport(2)
  }

  # Legends
  if (withlegend == "right" || withlegend == "bottom") {
    # Grid for legends
    if (withlegend == "right") {
      upViewport()
      downViewport("vplegend")
      pushViewport(viewport(layout = grid.layout(nrow = nplots, ncol = 1)))
      lposrow <- 1:nplots
      lposcol <- rep(1, nplots)
    } else {
      upViewport()
      downViewport("vplegend")
      pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = nplots)))
      lposrow <- rep(1, nplots)
      lposcol <- 1:nplots
    }

    # Legends for observations
    if (plots == "both" || plots == "obs") {
      for (i in 1:nchannels) {
        pushViewport(viewport(layout.pos.row = lposrow[i],
          layout.pos.col = lposcol[i]))
        pushViewport(viewport(width = unit(0.9, "npc")))
        par(plt = gridPLT(), new = TRUE)
        if (nchannels > 1) {
          if (withlegend == "bottom") {
            seqlegend(obs[[i]], fontsize = cex.legend, position = "top",
              ncol = ncol.legend[i], with.missing = with.missing.legend)
          } else {
            seqlegend(obs[[i]], fontsize = cex.legend, position = "left",
              ncol = ncol.legend[i], with.missing = with.missing.legend)
          }
        } else {
          if(withlegend == "bottom"){
            seqlegend(obs, fontsize = cex.legend, position = "top",
              ncol = ncol.legend[i], with.missing = with.missing.legend)
          } else {
            seqlegend(obs, fontsize = cex.legend, position = "left",
              ncol = ncol.legend[i], with.missing = with.missing.legend)
          }
        }
        popViewport(2)
      }
    }
    # Legends for most probable paths
    if (plots == "both" || plots == "hidden.paths") {
      pushViewport(viewport(layout.pos.row = lposrow[nplots],
        layout.pos.col = lposcol[nplots]))
      pushViewport(viewport(width = unit(0.9, "npc")))
      par(plt = gridPLT(), new = TRUE)
      if (withlegend == "bottom") {
        seqlegend(hidden.paths.seq, fontsize = cex.legend, position = "top",
          ncol = ncol.legend[length(ncol.legend)], with.missing = with.missing.legend)
      } else {
        seqlegend(hidden.paths.seq, fontsize = cex.legend, position = "left",
          ncol = ncol.legend[length(ncol.legend)],
          with.missing = with.missing.legend)
      }
      popViewport(2)
    }
    popViewport(2)
    # Combined legends
  } else if (withlegend == "right.combined" || withlegend == "bottom.combined") {
    ltext <- NULL
    cpal <- NULL
    if (plots == "both" || plots == "obs") {
      if (nchannels > 1) {
        for(i in 1:nchannels){
          ltext <- c(ltext, attr(obs[[i]], "labels"))
          cpal <- c(cpal, attr(obs[[i]], "cpal"))
        }
      } else {
        ltext <- c(ltext, attr(obs, "labels"))
        cpal <- c(cpal, attr(obs, "cpal"))
      }
    }
    if (plots == "both" || plots == "hidden.paths") {
      ltext <- c(ltext, attr(hidden.paths.seq, "labels"))
      cpal <- c(cpal, attr(hidden.paths.seq, "cpal"))
    }
    anymissing <- FALSE
    if (nchannels > 1) {
      for (i in 1:nchannels) {
        if (any(obs[[i]] == "*")) {
          anymissing <- TRUE
          break()
        }
      }
    } else {
      if (any(obs == "*")) {
        anymissing <- TRUE
      }
    }
    upViewport()
    downViewport("vplegend")
    par(plt = gridPLT(), new = TRUE)
    pushViewport(viewport(width = unit(0.9, "npc")))
    if (plots == "both" || plots == "obs") {
      if (nchannels > 1) {
        if (withlegend == "right.combined") {
          seqlegend(obs[[1]], fontsize = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = attr(obs[[1]],"missing.color"))
        } else { # withlegend == "bottom.combined"
          seqlegend(obs[[1]], fontsize = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = attr(obs[[1]],"missing.color"))
        }
      } else {
        if (withlegend == "right.combined") {
          seqlegend(obs, fontsize = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = attr(obs,"missing.color"))
        } else { # withlegend == "bottom.combined"
          seqlegend(obs, fontsize = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = attr(obs,"missing.color"))
        }
      }
    } else {
      if (nchannels > 1) {
        if (withlegend == "right.combined") {
          seqlegend(hidden.paths, fontsize = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = attr(obs[[1]],"missing.color"))
        } else { # withlegend == "bottom.combined"
          seqlegend(hidden.paths, fontsize = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = attr(obs[[1]],"missing.color"))
        }
      } else {
        if(withlegend == "right.combined"){
          seqlegend(hidden.paths, fontsize = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = attr(obs,"missing.color"))
        } else { # withlegend == "bottom.combined"
          seqlegend(hidden.paths, fontsize = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = attr(obs,"missing.color"))
        }
      }
    }
    popViewport(2)
  }

  # Title
  if (!is.na(title)) {
    if (title != FALSE) {
      upViewport()
      downViewport("vptitle")
      if (title.n) {
        if (nchannels > 1) {
          grid.text(paste0(title,", n = ", n.seq), y = unit(title.pos, "lines"),
            gp = gpar(cex = cex.title))
        } else {
          grid.text(paste0(title,", n = ", n.seq), y = unit(title.pos, "lines"),
            gp = gpar(cex = cex.title))
        }
        popViewport()
      } else {
        grid.text(title, y = unit(title.pos, "lines"),
          gp = gpar(cex = cex.title))
        popViewport()
      }
    }
  } else if (title.n) {
    upViewport()
    downViewport("vptitle")
    if (nchannels > 1) {
      grid.text(paste0("n = ", n.seq), y = unit(title.pos, "lines"),
        gp = gpar(cex = cex.title))
    } else {
      grid.text(paste0("n = ", n.seq), y = unit(title.pos, "lines"),
        gp = gpar(cex = cex.title))
    }
    popViewport()
  }
  popViewport()
}
