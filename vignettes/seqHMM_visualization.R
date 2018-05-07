## ----settingdata, message=FALSE, cache=FALSE, echo = FALSE, eval = TRUE----
library("seqHMM")

data("biofam", package = "TraMineR")
biofam_seq <- seqdef(biofam[, 10:25], start = 15, labels = c("parent",
  "left", "married", "left+marr", "child", "left+child", "left+marr+ch",
  "divorced"))

data("biofam3c")
marr_seq <- seqdef(biofam3c$married, start = 15, alphabet = c("single",
  "married", "divorced"))
child_seq <- seqdef(biofam3c$children, start = 15,
  alphabet = c("childless", "children"))
left_seq <- seqdef(biofam3c$left, start = 15, alphabet = c("with parents",
  "left home"))

attr(marr_seq, "cpal") <- c("violetred2", "darkgoldenrod2", "darkmagenta")
attr(child_seq, "cpal") <- c("darkseagreen1", "coral3")
attr(left_seq, "cpal") <- c("lightblue", "red3")

## ----seqplotSingle, fig.width=6.5, fig.height=3, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Sequence index plot (left) and state distribution plot (right) of annual family states for 100 individuals from the \\texttt{biofam} data.", fig.align='center', fig.keep='last', cache = TRUE, message=FALSE----

seqIplot <- ssp(
  biofam_seq[1:100, ], type = "I", sortv = "mds.obs", sort.channel = 1,
  xtlab = 15:30, xlab = "Age", 
  title.n = FALSE, ylab = "", title = "Sequences", title.pos = 0.5,
  legend.prop = 0.63, with.legend = FALSE)


seqdplot <- ssp(
  biofam_seq[1:100, ], type = "d",
  xtlab = 15:30, xlab = "Age", yaxis = TRUE, ylab = "Proportion",
  title.n = FALSE, title = "State distributions",
  legend.prop = 0.63, with.legend = FALSE)

gridplot(list(seqIplot, seqdplot), ncol = 3, col.prop = c(0.35, 0.4, 0.25))

## ----seqplotMulti, fig.width=6.5, fig.height=3.7, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Stacked sequence plot of the first ten individuals in the \\texttt{biofam} data. The top plot shows the original sequences, and the three bottom plots show the sequences in the separate channels for the same individuals. The sequences are in the same order in each plot, i.e., the same row always matches the same individual.", fig.align='center', fig.keep='last', cache = FALSE, message=FALSE----
seq_data <- list(
  biofam_seq[1:10,], marr_seq[1:10,], child_seq[1:10,],
  left_seq[1:10,])
ssplot(
  seq_data, type = "I", sortv = "from.start", sort.channel = 1,
  ylab = c("Original", "Marriage", "Parenthood", "Residence"),
  xtlab = 15:30, xlab = "Age", ylab.pos = c(1, 1.5), 
  title.n = FALSE, title = "Ten first sequences", 
  legend.prop = 0.63,
  ncol.legend = c(3, 1, 1, 1))

## ----HMMplot, fig.width=6.5, fig.height=4.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Illustrating a left-to-right hidden Markov model for the multichannel \\texttt{biofam} data as a directed graph. Pies represent the hidden states, with emission probabilities of combined observations as slices. Arrows illustrate transition probabilities between the hidden states. Probabilities of starting in each state are shown next to the pies.", fig.align='center', fig.keep='last', cache = TRUE, message=FALSE----
plot(hmm_biofam,
  layout = matrix(c(1, 2, 3, 4, 2,
                    1, 1, 1, 1, 0), ncol = 2),
  # varying curvature of edges
  edge.curved = c(0, -0.8, 0.6, 0, 0, -0.8, 0),
  # thinner edges and arrows
  cex.edge.width = 0.8, edge.arrow.size = 1,
  # fixing axes to the right scale
  xlim = c(0.5, 4.5), ylim = c(-0.5, 1.5), rescale = FALSE,
  # different legend properties
  with.legend = "bottom", legend.prop = 0.3, ncol.legend = 2,
  # distance of vertex labels to vertices
  vertex.label.dist = 1.1,
  # threshold for emission probabilities not shown as separate slices
  combine.slices = 0.02, combined.slice.label = "others (emission prob. < 0.02)"
  )

## ----settingsequences, message = FALSE-----------------------------------
library("seqHMM")

data("biofam3c")
marr_seq <- seqdef(
  biofam3c$married, start = 15, 
  alphabet = c("single", "married", "divorced"))
child_seq <- seqdef(
  biofam3c$children, start = 15, 
  alphabet = c("childless", "children"))
left_seq <- seqdef(
  biofam3c$left, start = 15, 
  alphabet = c("with parents", "left home"))

attr(marr_seq, "cpal") <- c("violetred2", "darkgoldenrod2", 
                            "darkmagenta")
attr(child_seq, "cpal") <- c("darkseagreen1", "coral3")
attr(left_seq, "cpal") <- c("lightblue", "red3")


## ----plottingsequences, fig.width=5, fig.height=3, dev.args=list(pointsize=10), fig.cap="Stacked sequence plot of annual state distributions in the three-channel \\texttt{biofam} data. This is the default output of the \\texttt{ssplot} function.", fig.keep='last', fig.align='center', cache=FALSE, echo = TRUE----

ssplot(
  x = list("Marriage" = marr_seq, "Parenthood" = child_seq,
       "Residence" = left_seq))


## ----plottingsequencesHMM, fig.width=5, fig.height=5, dev.args=list(pointsize=10), fig.cap="Stacked sequence plot of observations and hidden state paths using a hidden Markov model object.", fig.keep='last', fig.align='center', cache=FALSE, echo = TRUE----

data("hmm_biofam")

ssplot(x = hmm_biofam, plots = "both")


## ----seqIplot, fig.width=5, fig.height=3, dev.args=list(pointsize=10), fig.cap="Sequence index plot showing observed sequences sorted by the third channel, residence.", fig.keep='last', fig.align='center', cache=TRUE, echo = TRUE----

ssplot(
  hmm_biofam, type = "I", sortv = "from.start", sort.channel = 3)


## ----ssp, fig.width=6, fig.height=3.2, dev.args=list(pointsize=10), fig.cap="Example on saving \\texttt{ssp} objects. Sequences are sorted according to multidimensional scaling scores.", fig.keep='last', fig.align='center', cache=TRUE, echo = TRUE----

ssp_def <- ssp(
  hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden", 
  ylab.pos = c(1, 2), 
  title = "Family trajectories", title.n = FALSE,
  xtlab = 15:30, xlab = "Age",
  ncol.legend = c(2, 1, 1), legend.prop = 0.37)

plot(ssp_def)


## ----gridplot1, fig.width=5.5, fig.height=3.5, dev.args=list(pointsize=10), echo=TRUE, fig.cap="Showing state distribution plots for women and men in the \\texttt{biofam} data. Two figures were defined with the \\texttt{ssp} function and then combined into one figure with the \\texttt{gridplot} function.", fig.align='center', fig.keep='last', cache = TRUE----
ssp_f <- ssp(
  list(marr_seq[biofam3c$covariates$sex == "woman",],
       child_seq[biofam3c$covariates$sex == "woman",],
       left_seq[biofam3c$covariates$sex == "woman",]),
  type = "I", sortv = "mds.obs", with.legend = FALSE, 
  title = "Women", xtlab = 15:30, ylab.pos = c(1, 2, 1), 
  ylab = c("Married", "Children", "Residence"))

ssp_m <- update(
  ssp_f, title = "Men",
  x = list(marr_seq[biofam3c$covariates$sex == "man",],
           child_seq[biofam3c$covariates$sex == "man",],
           left_seq[biofam3c$covariates$sex == "man",]))

gridplot(
  list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
  legend.pos = "bottom", legend.pos2 = "top", 
  row.prop = c(0.65, 0.35))

## ----gridplot2, fig.width=5.5, fig.height=4.5, dev.args=list(pointsize=10), echo=TRUE, fig.cap="Another example of \\texttt{gridplot}. Showing sequences and state distributions for women and men.", fig.align='center', fig.keep='last', cache = TRUE----

ssp_f2 <- update(ssp_f, type = "d", title = FALSE)

ssp_m2 <- update(ssp_m, type = "d", title = FALSE)

gridplot(
  list(ssp_f, ssp_f2, ssp_m, ssp_m2), ncol = 3, nrow = 2,
  legend.pos = 3:4)

## ----gridplot3, fig.width=5.5, fig.height=4, dev.args=list(pointsize=10), echo=TRUE, fig.cap="Another example of \\texttt{gridplot}. Showing three-channel sequences and state distributions of combined states for women.", fig.align='center', fig.keep='last', cache = TRUE----

ssp_f3 <- update(
  ssp_f, with.legend = TRUE, legend.prop = 0.4, ylab.pos = 1, 
  cex.lab = 0.9, cex.axis = 0.8, cex.legend = 0.9)

ssp_f4 <- ssp(
  biofam_seq[biofam3c$covariates$sex == "woman",],
  type = "d", title.n = FALSE, xtlab = 15:30, 
  title = "State distributions for \n combined states (women)", 
  title.pos = 1.5, ylab = "", xlab = "Age",
  with.legend = "bottom", ncol.legend = 2, 
  cex.lab = 0.9, cex.axis = 0.8, cex.legend = 0.9)

gridplot(list(ssp_f3, ssp_f4), with.legend = FALSE, ncol = 2, 
         col.prop = c(0.55, 0.45))

## ----code_plottingHMMbasic, fig.height=5, fig.width=8, echo=TRUE, fig.align='center', fig.keep='last', cache = TRUE, eval = TRUE, fig.cap="A default plot of a hidden Markov model."----
plot(hmm_biofam)

## ----HMMplotCode, echo=TRUE, eval=FALSE----------------------------------
#  plot(hmm_biofam,
#    layout = matrix(c(1, 2, 3, 4, 2,
#                      1, 1, 1, 1, 0), ncol = 2),
#    xlim = c(0.5, 4.5), ylim = c(-0.5, 1.5), rescale = FALSE,
#    edge.curved = c(0, -0.8, 0.6, 0, 0, -0.8, 0),
#    cex.edge.width = 0.8, edge.arrow.size = 1,
#    legend.prop = 0.3, ncol.legend = 2,
#    vertex.label.dist = 1.1, combine.slices = 0.02,
#    combined.slice.label = "others (emission prob. < 0.02)")

## ----HMMplotLayout, fig.width=5.5, fig.height=4, dev.args=list(pointsize=10), echo=TRUE, fig.cap="Another example of \\texttt{plot.hmm}.", fig.align='center', fig.keep='last', cache = TRUE, message=FALSE----

require("igraph")
set.seed(1234)
plot(hmm_biofam,
  layout = layout_nicely, pie = FALSE, 
  vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
  edge.curved = FALSE, edge.width = 1, 
  loops = TRUE, edge.loop.angle = -pi/8,
  trim = 0.01, label.signif = 3,
  xlim = c(-1, 1.3))

## ----colorpalette, fig.width=5.5, fig.height=3, dev.args=list(pointsize=10), echo=TRUE, fig.cap="Helper function for plotting colour palettes with their names.", fig.align='center', fig.keep='last', cache = TRUE----

plot_colors(colorpalette[[7]])


