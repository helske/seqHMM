## ----setup, include=FALSE, cache=FALSE------------------------------------
library(knitr)
opts_chunk$set(concordance = TRUE, tidy = FALSE)
options(prompt = "R> ", continue = "+  ", width = 76, useFancyQuotes = FALSE)


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

## ----graphicalillustrations2, fig.width=6.5, fig.height=3.7, dev.args=list(pointsize=10), fig.keep='last', cache=FALSE, message=FALSE, echo=FALSE, fig.cap='Stacked sequence plot of the first ten individuals in the \\code{biofam} data plotted with the \\code{ssplot} function. The top plot shows the original sequences, and the three bottom plots show the sequences in the separate channels for the same individuals. The sequences are in the same order in each plot, i.e., the same row always matches the same individual.', fig.align='center'----
ssplot(list(biofam_seq[1:10,], marr_seq[1:10,], child_seq[1:10,],
  left_seq[1:10,]),
  sortv = "from.start", sort.channel = 1, type = "I",
  ylab = c("Original", "Marriage", "Parenthood", "Residence"),
  xtlab = 15:30, xlab = "Age", title = "Ten first sequences",
  title.n = FALSE, legend.prop = 0.63, ylab.pos = c(1, 1.5),
  ncol.legend = c(3, 1, 1, 1))

## ----plottingsequences, fig.width=5, fig.height=3, dev.args=list(pointsize=10), fig.cap="Stacked sequence plot of annual state distributions in the three-channel \\code{biofam} data. This is the default output of the \\code{ssplot} function. The labels for the channels are taken from the named list of state sequence objects, and the labels for the x axis ticks are taken from the column names of the first object.", fig.keep='last', fig.align='center', cache=FALSE, echo = FALSE----
ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq,
  "Residence" = left_seq))

## ----gridplot1, fig.width=5.5, fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Showing state distribution plots for women and men in the \\code{biofam} data. Two figures were defined with the \\code{ssp} function and then combined into one figure with the \\code{gridplot} function.", fig.align='center', fig.keep='last', cache = FALSE----
ssp_f <- ssp(list(marr_seq[biofam3c$covariates$sex == "woman",],
    child_seq[biofam3c$covariates$sex == "woman",],
    left_seq[biofam3c$covariates$sex == "woman",]),
  type = "I", sortv = "mds.obs", with.legend = FALSE, title = "Women", 
  ylab.pos = c(1, 2, 1), xtlab = 15:30, ylab = c("Married", "Children", 
    "Residence"))

ssp_m <- update(ssp_f, title = "Men", 
  x = list(marr_seq[biofam3c$covariates$sex == "man",],
    child_seq[biofam3c$covariates$sex == "man",],
    left_seq[biofam3c$covariates$sex == "man",]))

gridplot(list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
  legend.pos = "bottom", legend.pos2 = "top", row.prop = c(0.65, 0.35))

## ----code_mcHMM, cache=FALSE, echo = FALSE, message=FALSE, warning=TRUE, eval = TRUE----
mc_init <- c(0.9, 0.05, 0.02, 0.02, 0.01)

mc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0, 0.90, 0.05, 
  0.03, 0.02, 0, 0, 0.90, 0.07, 0.03, 0, 0, 0, 0.90, 0.10, 0, 0, 0,    
  0, 1), nrow = 5, ncol = 5, byrow = TRUE)

mc_emiss_marr <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 
  0.90, 0.05, 0.05, 0.90, 0.05, 0.30, 0.30, 0.40), nrow = 5, ncol = 3, 
  byrow = TRUE)

mc_emiss_child <- matrix(c(0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9,
  0.5, 0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_emiss_left <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9,
  0.5, 0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_obs <- list(marr_seq, child_seq, left_seq)

mc_emiss <- list(mc_emiss_marr, mc_emiss_child, mc_emiss_left)

mc_initmod <- build_hmm(observations = mc_obs, initial_probs = mc_init, 
  transition_probs = mc_trans, emission_probs = mc_emiss,
  channel_names = c("Marriage", "Parenthood", "Residence"))

# For CRAN vignette: load the estimated model object for speed-up
data("hmm_biofam")
# mc_fit <- fit_model(mc_initmod, em_step = FALSE, local_step = TRUE,
#   threads = 4)


## ----plottingHMM, out.width='\\linewidth', fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Illustrating a hidden Markov model as a directed graph. Pies represent five hidden states, with slices showing emission probabilities of combinations of observed states. States with emission probability less than 0.05 are combined into one slice. Edges show the transtion probabilities. Initial probabilities of hidden states are given below the pies.", fig.align='center', fig.keep='last', cache = FALSE----
plot(hmm_biofam, vertex.size = 50, vertex.label.dist = 1.5,
  edge.curved = c(0, 0.6, -0.8, 0.6, 0, 0.6, 0), legend.prop = 0.3, 
  combined.slice.label = "States with prob. < 0.05")

## ----graphicalillustrations5, out.width='\\linewidth', fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Another version of the hidden Markov model of Figure 4 with a different layout and modified labels, legends, and colors. All observed states are shown.", fig.align='center', fig.keep='last', cache = FALSE----
vertex_layout <- matrix(c(1, 2, 2, 3, 1, 0, 0.5, -0.5, 0, -1), 
  ncol = 2)

plot(hmm_biofam, layout = vertex_layout, xlim = c(0.5, 3.5), 
  ylim = c(-1.5, 1), rescale = FALSE, vertex.size = 50, 
  vertex.label.pos = c("left", "top", "bottom", "right", "left"),
  edge.curved = FALSE, edge.width = 1, edge.arrow.size = 1, 
  with.legend = "left", legend.prop = 0.4, label.signif = 1, 
  combine.slices = 0, cpal = colorpalette[[30]][c(14:5)])

## ----ssplotHMM, fig.width=5.5, fig.height=5.5, dev.args=list(pointsize=10), fig.cap="Using the \\code{ssplot} function for an \\code{hmm} object makes it easy to plot the observed sequences together with the most probable paths of hidden states given the model.", fig.align='center', fig.keep='last', cache = FALSE, echo = FALSE----
ssplot(hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden",
  title = "Observed and hidden state sequences", xtlab = 15:30, 
  xlab = "Age")

## ----code_settingdata, ref.label = 'settingdata', message=FALSE, warning = TRUE, echo = TRUE----
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

## ----code_plottingsequences, ref.label = 'plottingsequences', echo = TRUE, eval = FALSE----
#  ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq,
#    "Residence" = left_seq))

## ----code_graphicalillustrations2, ref.label = 'graphicalillustrations2', echo=TRUE, eval = FALSE----
#  ssplot(list(biofam_seq[1:10,], marr_seq[1:10,], child_seq[1:10,],
#    left_seq[1:10,]),
#    sortv = "from.start", sort.channel = 1, type = "I",
#    ylab = c("Original", "Marriage", "Parenthood", "Residence"),
#    xtlab = 15:30, xlab = "Age", title = "Ten first sequences",
#    title.n = FALSE, legend.prop = 0.63, ylab.pos = c(1, 1.5),
#    ncol.legend = c(3, 1, 1, 1))

## ----code_gridplot1, ref.label = 'gridplot1', echo=TRUE, eval = FALSE-----
#  ssp_f <- ssp(list(marr_seq[biofam3c$covariates$sex == "woman",],
#      child_seq[biofam3c$covariates$sex == "woman",],
#      left_seq[biofam3c$covariates$sex == "woman",]),
#    type = "I", sortv = "mds.obs", with.legend = FALSE, title = "Women",
#    ylab.pos = c(1, 2, 1), xtlab = 15:30, ylab = c("Married", "Children",
#      "Residence"))
#  
#  ssp_m <- update(ssp_f, title = "Men",
#    x = list(marr_seq[biofam3c$covariates$sex == "man",],
#      child_seq[biofam3c$covariates$sex == "man",],
#      left_seq[biofam3c$covariates$sex == "man",]))
#  
#  gridplot(list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
#    legend.pos = "bottom", legend.pos2 = "top", row.prop = c(0.65, 0.35))

## ----code_sc_buildHMM_random, cache=FALSE---------------------------------
sc_initmod_random <- build_hmm(observations = biofam_seq, n_states = 5)

## ----code_sc_initialvalues, cache=FALSE-----------------------------------
sc_init <- c(0.9, 0.06, 0.02, 0.01, 0.01)

sc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0.02, 0.80, 0.10, 
  0.05, 0.03, 0.02, 0.03, 0.80, 0.10, 0.05, 0.02, 0.03, 0.05, 0.80, 0.10, 
  0.02, 0.03, 0.05, 0.05, 0.85), nrow = 5, ncol = 5, byrow = TRUE)

sc_emiss <- matrix(NA, nrow = 5, ncol = 8)
sc_emiss[1,] <- seqstatf(biofam_seq[, 1:4])[, 2] + 0.1
sc_emiss[2,] <- seqstatf(biofam_seq[, 5:7])[, 2] + 0.1
sc_emiss[3,] <- seqstatf(biofam_seq[, 8:10])[, 2] + 0.1
sc_emiss[4,] <- seqstatf(biofam_seq[, 11:13])[, 2] + 0.1
sc_emiss[5,] <- seqstatf(biofam_seq[, 14:16])[, 2] + 0.1
sc_emiss <- sc_emiss / rowSums(sc_emiss)

rownames(sc_trans) <- colnames(sc_trans) <- rownames(sc_emiss) <-
  paste("State", 1:5)

colnames(sc_emiss) <- attr(biofam_seq, "labels")

sc_trans
round(sc_emiss, 3)

## ----code_sc_buildHMM, cache=FALSE----------------------------------------
sc_initmod <- build_hmm(observations = biofam_seq, initial_probs = sc_init,
  transition_probs = sc_trans, emission_probs = sc_emiss)

## ----code_sc_fitHMM, cache=FALSE------------------------------------------
sc_fit <- fit_model(sc_initmod)

## ----code_sc_results1, cache=FALSE----------------------------------------
sc_fit$logLik

## ----code_sc_results2, cache=FALSE----------------------------------------
sc_fit$model

## ----code_mcHMM2, echo = TRUE, message=TRUE, warnings=TRUE, eval = TRUE, cache = FALSE----
mc_init <- c(0.9, 0.05, 0.02, 0.02, 0.01)

mc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0, 0.90, 0.05, 0.03, 
  0.02, 0, 0, 0.90, 0.07, 0.03, 0, 0, 0, 0.90, 0.10, 0, 0, 0, 0, 1), 
  nrow = 5, ncol = 5, byrow = TRUE)

mc_emiss_marr <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 0.90, 
  0.05, 0.05, 0.90, 0.05, 0.30, 0.30, 0.40), nrow = 5, ncol = 3, 
  byrow = TRUE)

mc_emiss_child <- matrix(c(0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.5, 
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_emiss_left <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.5, 
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_obs <- list(marr_seq, child_seq, left_seq)

mc_emiss <- list(mc_emiss_marr, mc_emiss_child, mc_emiss_left)

mc_initmod <- build_hmm(observations = mc_obs, initial_probs = mc_init, 
  transition_probs = mc_trans, emission_probs = mc_emiss,
  channel_names = c("Marriage", "Parenthood", "Residence"))

# For CRAN vignette: load the estimated model object for speed-up
data("hmm_biofam")
# mc_fit <- fit_model(mc_initmod, em_step = FALSE, local_step = TRUE,
# threads = 4)


## ----code_mcHMM_BIC, cache=FALSE, echo = TRUE, message=FALSE, eval = TRUE----
# Vignette: already loaded hmm_biofam
# hmm_biofam <- mc_fit$model
BIC(hmm_biofam)

## ----code_MHMM, cache=FALSE, echo = TRUE, eval = TRUE, warning=TRUE-------
mc_init2 <- c(0.9, 0.05, 0.03, 0.02)

mc_trans2 <- matrix(c(0.85, 0.05, 0.05, 0.05, 0, 0.90, 0.05, 0.05, 0, 0, 
  0.95, 0.05, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

mc_emiss_marr2 <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 
  0.85, 0.10, 0.05, 0.80, 0.15), nrow = 4, ncol = 3, byrow = TRUE)

mc_emiss_child2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mc_emiss_left2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mhmm_init <- list(mc_init, mc_init2)

mhmm_trans <- list(mc_trans, mc_trans2)

mhmm_emiss <- list(list(mc_emiss_marr, mc_emiss_child, mc_emiss_left), 
  list(mc_emiss_marr2, mc_emiss_child2, mc_emiss_left2))

biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, 
  c(1908, 1935, 1945, 1957))
biofam3c$covariates$cohort <- factor(biofam3c$covariates$cohort, 
  labels=c("1909-1935", "1936-1945", "1946-1957"))

init_mhmm <- build_mhmm(observations = mc_obs, initial_probs = mhmm_init, 
  transition_probs = mhmm_trans, emission_probs = mhmm_emiss, 
  formula = ~sex + cohort, data = biofam3c$covariates, 
  channel_names = c("Marriage", "Parenthood", "Residence"), 
  cluster_names = c("Cluster 1", "Cluster 2"))

# vignette: less restarts and no parallelization
set.seed(1011)
mhmm_fit <- fit_model(init_mhmm, local_step = TRUE, threads = 1,
  control_em = list(restart = list(times = 10)))
mhmm <- mhmm_fit$model

## ----code_summaryMHMM, cache=FALSE----------------------------------------
summary(mhmm, conditional_se = FALSE)

## ----code_plottingHMMbasic, out.width='\\linewidth', fig.height=4, dev.args=list(pointsize=10), echo=TRUE, fig.align='center', fig.keep='last', cache = FALSE, eval = TRUE, fig.cap="A default plot of a hidden Markov model."----
plot(hmm_biofam)

## ----code_plottingHMM, ref.label='plottingHMM', echo=TRUE, eval = FALSE----
#  plot(hmm_biofam, vertex.size = 50, vertex.label.dist = 1.5,
#    edge.curved = c(0, 0.6, -0.8, 0.6, 0, 0.6, 0), legend.prop = 0.3,
#    combined.slice.label = "States with prob. < 0.05")

## ----code_graphicalillustrations5, ref.label = 'graphicalillustrations5', echo=TRUE, eval = FALSE----
#  vertex_layout <- matrix(c(1, 2, 2, 3, 1, 0, 0.5, -0.5, 0, -1),
#    ncol = 2)
#  
#  plot(hmm_biofam, layout = vertex_layout, xlim = c(0.5, 3.5),
#    ylim = c(-1.5, 1), rescale = FALSE, vertex.size = 50,
#    vertex.label.pos = c("left", "top", "bottom", "right", "left"),
#    edge.curved = FALSE, edge.width = 1, edge.arrow.size = 1,
#    with.legend = "left", legend.prop = 0.4, label.signif = 1,
#    combine.slices = 0, cpal = colorpalette[[30]][c(14:5)])

## ----code_ssplotHMM, ref.label = 'ssplotHMM', eval = FALSE, echo = TRUE----
#  ssplot(hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden",
#    title = "Observed and hidden state sequences", xtlab = 15:30,
#    xlab = "Age")

## ----code_plottingMHMMbasic, fig.width=6.5, fig.height=8, dev.args=list(pointsize=10), echo=TRUE, fig.align='center', fig.keep='last', cache = FALSE, eval = TRUE, fig.cap="Plotting submodels of an MHMM with the \\code{plot} method."----
plot(mhmm, interactive = FALSE, nrow = 2, legend.prop = 0.45,
  vertex.size = 50, vertex.label.cex = 1.3, cex.legend = 1.3, 
  edge.curved = 0.65, edge.label.cex = 1.3, edge.arrow.size = 0.8)

## ----code_ssplotMHMM, eval = FALSE, echo = TRUE---------------------------
#  mssplot(mhmm, ask = TRUE)

## ----setup, include=FALSE, cache=FALSE------------------------------------
library(knitr)
opts_chunk$set(concordance = TRUE, tidy = FALSE)
options(prompt = "R> ", continue = "+  ", width = 76, useFancyQuotes = FALSE)


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

## ----graphicalillustrations2, fig.width=6.5, fig.height=3.7, dev.args=list(pointsize=10), fig.keep='last', cache=FALSE, message=FALSE, echo=FALSE, fig.cap='Stacked sequence plot of the first ten individuals in the \\code{biofam} data plotted with the \\code{ssplot} function. The top plot shows the original sequences, and the three bottom plots show the sequences in the separate channels for the same individuals. The sequences are in the same order in each plot, i.e., the same row always matches the same individual.', fig.align='center'----
ssplot(list(biofam_seq[1:10,], marr_seq[1:10,], child_seq[1:10,],
  left_seq[1:10,]),
  sortv = "from.start", sort.channel = 1, type = "I",
  ylab = c("Original", "Marriage", "Parenthood", "Residence"),
  xtlab = 15:30, xlab = "Age", title = "Ten first sequences",
  title.n = FALSE, legend.prop = 0.63, ylab.pos = c(1, 1.5),
  ncol.legend = c(3, 1, 1, 1))

## ----plottingsequences, fig.width=5, fig.height=3, dev.args=list(pointsize=10), fig.cap="Stacked sequence plot of annual state distributions in the three-channel \\code{biofam} data. This is the default output of the \\code{ssplot} function. The labels for the channels are taken from the named list of state sequence objects, and the labels for the x axis ticks are taken from the column names of the first object.", fig.keep='last', fig.align='center', cache=FALSE, echo = FALSE----
ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq,
  "Residence" = left_seq))

## ----gridplot1, fig.width=5.5, fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Showing state distribution plots for women and men in the \\code{biofam} data. Two figures were defined with the \\code{ssp} function and then combined into one figure with the \\code{gridplot} function.", fig.align='center', fig.keep='last', cache = FALSE----
ssp_f <- ssp(list(marr_seq[biofam3c$covariates$sex == "woman",],
    child_seq[biofam3c$covariates$sex == "woman",],
    left_seq[biofam3c$covariates$sex == "woman",]),
  type = "I", sortv = "mds.obs", with.legend = FALSE, title = "Women", 
  ylab.pos = c(1, 2, 1), xtlab = 15:30, ylab = c("Married", "Children", 
    "Residence"))

ssp_m <- update(ssp_f, title = "Men", 
  x = list(marr_seq[biofam3c$covariates$sex == "man",],
    child_seq[biofam3c$covariates$sex == "man",],
    left_seq[biofam3c$covariates$sex == "man",]))

gridplot(list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
  legend.pos = "bottom", legend.pos2 = "top", row.prop = c(0.65, 0.35))

## ----code_mcHMM, cache=FALSE, echo = FALSE, message=FALSE, warning=TRUE, eval = TRUE----
mc_init <- c(0.9, 0.05, 0.02, 0.02, 0.01)

mc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0, 0.90, 0.05, 
  0.03, 0.02, 0, 0, 0.90, 0.07, 0.03, 0, 0, 0, 0.90, 0.10, 0, 0, 0,    
  0, 1), nrow = 5, ncol = 5, byrow = TRUE)

mc_emiss_marr <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 
  0.90, 0.05, 0.05, 0.90, 0.05, 0.30, 0.30, 0.40), nrow = 5, ncol = 3, 
  byrow = TRUE)

mc_emiss_child <- matrix(c(0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9,
  0.5, 0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_emiss_left <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9,
  0.5, 0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_obs <- list(marr_seq, child_seq, left_seq)

mc_emiss <- list(mc_emiss_marr, mc_emiss_child, mc_emiss_left)

mc_initmod <- build_hmm(observations = mc_obs, initial_probs = mc_init, 
  transition_probs = mc_trans, emission_probs = mc_emiss,
  channel_names = c("Marriage", "Parenthood", "Residence"))

# For CRAN vignette: load the estimated model object for speed-up
data("hmm_biofam")
# mc_fit <- fit_model(mc_initmod, em_step = FALSE, local_step = TRUE,
#   threads = 4)


## ----plottingHMM, out.width='\\linewidth', fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Illustrating a hidden Markov model as a directed graph. Pies represent five hidden states, with slices showing emission probabilities of combinations of observed states. States with emission probability less than 0.05 are combined into one slice. Edges show the transtion probabilities. Initial probabilities of hidden states are given below the pies.", fig.align='center', fig.keep='last', cache = FALSE----
plot(hmm_biofam, vertex.size = 50, vertex.label.dist = 1.5,
  edge.curved = c(0, 0.6, -0.8, 0.6, 0, 0.6, 0), legend.prop = 0.3, 
  combined.slice.label = "States with prob. < 0.05")

## ----graphicalillustrations5, out.width='\\linewidth', fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Another version of the hidden Markov model of Figure 4 with a different layout and modified labels, legends, and colors. All observed states are shown.", fig.align='center', fig.keep='last', cache = FALSE----
vertex_layout <- matrix(c(1, 2, 2, 3, 1, 0, 0.5, -0.5, 0, -1), 
  ncol = 2)

plot(hmm_biofam, layout = vertex_layout, xlim = c(0.5, 3.5), 
  ylim = c(-1.5, 1), rescale = FALSE, vertex.size = 50, 
  vertex.label.pos = c("left", "top", "bottom", "right", "left"),
  edge.curved = FALSE, edge.width = 1, edge.arrow.size = 1, 
  with.legend = "left", legend.prop = 0.4, label.signif = 1, 
  combine.slices = 0, cpal = colorpalette[[30]][c(14:5)])

## ----ssplotHMM, fig.width=5.5, fig.height=5.5, dev.args=list(pointsize=10), fig.cap="Using the \\code{ssplot} function for an \\code{hmm} object makes it easy to plot the observed sequences together with the most probable paths of hidden states given the model.", fig.align='center', fig.keep='last', cache = FALSE, echo = FALSE----
ssplot(hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden",
  title = "Observed and hidden state sequences", xtlab = 15:30, 
  xlab = "Age")

## ----code_settingdata, ref.label = 'settingdata', message=FALSE, warning = TRUE, echo = TRUE----
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

## ----code_plottingsequences, ref.label = 'plottingsequences', echo = TRUE, eval = FALSE----
#  ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq,
#    "Residence" = left_seq))

## ----code_graphicalillustrations2, ref.label = 'graphicalillustrations2', echo=TRUE, eval = FALSE----
#  ssplot(list(biofam_seq[1:10,], marr_seq[1:10,], child_seq[1:10,],
#    left_seq[1:10,]),
#    sortv = "from.start", sort.channel = 1, type = "I",
#    ylab = c("Original", "Marriage", "Parenthood", "Residence"),
#    xtlab = 15:30, xlab = "Age", title = "Ten first sequences",
#    title.n = FALSE, legend.prop = 0.63, ylab.pos = c(1, 1.5),
#    ncol.legend = c(3, 1, 1, 1))

## ----code_gridplot1, ref.label = 'gridplot1', echo=TRUE, eval = FALSE-----
#  ssp_f <- ssp(list(marr_seq[biofam3c$covariates$sex == "woman",],
#      child_seq[biofam3c$covariates$sex == "woman",],
#      left_seq[biofam3c$covariates$sex == "woman",]),
#    type = "I", sortv = "mds.obs", with.legend = FALSE, title = "Women",
#    ylab.pos = c(1, 2, 1), xtlab = 15:30, ylab = c("Married", "Children",
#      "Residence"))
#  
#  ssp_m <- update(ssp_f, title = "Men",
#    x = list(marr_seq[biofam3c$covariates$sex == "man",],
#      child_seq[biofam3c$covariates$sex == "man",],
#      left_seq[biofam3c$covariates$sex == "man",]))
#  
#  gridplot(list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
#    legend.pos = "bottom", legend.pos2 = "top", row.prop = c(0.65, 0.35))

## ----code_sc_buildHMM_random, cache=FALSE---------------------------------
sc_initmod_random <- build_hmm(observations = biofam_seq, n_states = 5)

## ----code_sc_initialvalues, cache=FALSE-----------------------------------
sc_init <- c(0.9, 0.06, 0.02, 0.01, 0.01)

sc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0.02, 0.80, 0.10, 
  0.05, 0.03, 0.02, 0.03, 0.80, 0.10, 0.05, 0.02, 0.03, 0.05, 0.80, 0.10, 
  0.02, 0.03, 0.05, 0.05, 0.85), nrow = 5, ncol = 5, byrow = TRUE)

sc_emiss <- matrix(NA, nrow = 5, ncol = 8)
sc_emiss[1,] <- seqstatf(biofam_seq[, 1:4])[, 2] + 0.1
sc_emiss[2,] <- seqstatf(biofam_seq[, 5:7])[, 2] + 0.1
sc_emiss[3,] <- seqstatf(biofam_seq[, 8:10])[, 2] + 0.1
sc_emiss[4,] <- seqstatf(biofam_seq[, 11:13])[, 2] + 0.1
sc_emiss[5,] <- seqstatf(biofam_seq[, 14:16])[, 2] + 0.1
sc_emiss <- sc_emiss / rowSums(sc_emiss)

rownames(sc_trans) <- colnames(sc_trans) <- rownames(sc_emiss) <-
  paste("State", 1:5)

colnames(sc_emiss) <- attr(biofam_seq, "labels")

sc_trans
round(sc_emiss, 3)

## ----code_sc_buildHMM, cache=FALSE----------------------------------------
sc_initmod <- build_hmm(observations = biofam_seq, initial_probs = sc_init,
  transition_probs = sc_trans, emission_probs = sc_emiss)

## ----code_sc_fitHMM, cache=FALSE------------------------------------------
sc_fit <- fit_model(sc_initmod)

## ----code_sc_results1, cache=FALSE----------------------------------------
sc_fit$logLik

## ----code_sc_results2, cache=FALSE----------------------------------------
sc_fit$model

## ----code_mcHMM2, echo = TRUE, message=TRUE, warnings=TRUE, eval = TRUE, cache = FALSE----
mc_init <- c(0.9, 0.05, 0.02, 0.02, 0.01)

mc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0, 0.90, 0.05, 0.03, 
  0.02, 0, 0, 0.90, 0.07, 0.03, 0, 0, 0, 0.90, 0.10, 0, 0, 0, 0, 1), 
  nrow = 5, ncol = 5, byrow = TRUE)

mc_emiss_marr <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 0.90, 
  0.05, 0.05, 0.90, 0.05, 0.30, 0.30, 0.40), nrow = 5, ncol = 3, 
  byrow = TRUE)

mc_emiss_child <- matrix(c(0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.5, 
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_emiss_left <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.5, 
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_obs <- list(marr_seq, child_seq, left_seq)

mc_emiss <- list(mc_emiss_marr, mc_emiss_child, mc_emiss_left)

mc_initmod <- build_hmm(observations = mc_obs, initial_probs = mc_init, 
  transition_probs = mc_trans, emission_probs = mc_emiss,
  channel_names = c("Marriage", "Parenthood", "Residence"))

# For CRAN vignette: load the estimated model object for speed-up
data("hmm_biofam")
# mc_fit <- fit_model(mc_initmod, em_step = FALSE, local_step = TRUE,
# threads = 4)


## ----code_mcHMM_BIC, cache=FALSE, echo = TRUE, message=FALSE, eval = TRUE----
# Vignette: already loaded hmm_biofam
# hmm_biofam <- mc_fit$model
BIC(hmm_biofam)

## ----code_MHMM, cache=FALSE, echo = TRUE, eval = TRUE, warning=TRUE-------
mc_init2 <- c(0.9, 0.05, 0.03, 0.02)

mc_trans2 <- matrix(c(0.85, 0.05, 0.05, 0.05, 0, 0.90, 0.05, 0.05, 0, 0, 
  0.95, 0.05, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

mc_emiss_marr2 <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 
  0.85, 0.10, 0.05, 0.80, 0.15), nrow = 4, ncol = 3, byrow = TRUE)

mc_emiss_child2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mc_emiss_left2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mhmm_init <- list(mc_init, mc_init2)

mhmm_trans <- list(mc_trans, mc_trans2)

mhmm_emiss <- list(list(mc_emiss_marr, mc_emiss_child, mc_emiss_left), 
  list(mc_emiss_marr2, mc_emiss_child2, mc_emiss_left2))

biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, 
  c(1908, 1935, 1945, 1957))
biofam3c$covariates$cohort <- factor(biofam3c$covariates$cohort, 
  labels=c("1909-1935", "1936-1945", "1946-1957"))

init_mhmm <- build_mhmm(observations = mc_obs, initial_probs = mhmm_init, 
  transition_probs = mhmm_trans, emission_probs = mhmm_emiss, 
  formula = ~sex + cohort, data = biofam3c$covariates, 
  channel_names = c("Marriage", "Parenthood", "Residence"), 
  cluster_names = c("Cluster 1", "Cluster 2"))

# vignette: less restarts and no parallelization
set.seed(1011)
mhmm_fit <- fit_model(init_mhmm, local_step = TRUE, threads = 1,
  control_em = list(restart = list(times = 10)))
mhmm <- mhmm_fit$model

## ----code_summaryMHMM, cache=FALSE----------------------------------------
summary(mhmm, conditional_se = FALSE)

## ----code_plottingHMMbasic, out.width='\\linewidth', fig.height=4, dev.args=list(pointsize=10), echo=TRUE, fig.align='center', fig.keep='last', cache = FALSE, eval = TRUE, fig.cap="A default plot of a hidden Markov model."----
plot(hmm_biofam)

## ----code_plottingHMM, ref.label='plottingHMM', echo=TRUE, eval = FALSE----
#  plot(hmm_biofam, vertex.size = 50, vertex.label.dist = 1.5,
#    edge.curved = c(0, 0.6, -0.8, 0.6, 0, 0.6, 0), legend.prop = 0.3,
#    combined.slice.label = "States with prob. < 0.05")

## ----code_graphicalillustrations5, ref.label = 'graphicalillustrations5', echo=TRUE, eval = FALSE----
#  vertex_layout <- matrix(c(1, 2, 2, 3, 1, 0, 0.5, -0.5, 0, -1),
#    ncol = 2)
#  
#  plot(hmm_biofam, layout = vertex_layout, xlim = c(0.5, 3.5),
#    ylim = c(-1.5, 1), rescale = FALSE, vertex.size = 50,
#    vertex.label.pos = c("left", "top", "bottom", "right", "left"),
#    edge.curved = FALSE, edge.width = 1, edge.arrow.size = 1,
#    with.legend = "left", legend.prop = 0.4, label.signif = 1,
#    combine.slices = 0, cpal = colorpalette[[30]][c(14:5)])

## ----code_ssplotHMM, ref.label = 'ssplotHMM', eval = FALSE, echo = TRUE----
#  ssplot(hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden",
#    title = "Observed and hidden state sequences", xtlab = 15:30,
#    xlab = "Age")

## ----code_plottingMHMMbasic, fig.width=6.5, fig.height=8, dev.args=list(pointsize=10), echo=TRUE, fig.align='center', fig.keep='last', cache = FALSE, eval = TRUE, fig.cap="Plotting submodels of an MHMM with the \\code{plot} method."----
plot(mhmm, interactive = FALSE, nrow = 2, legend.prop = 0.45,
  vertex.size = 50, vertex.label.cex = 1.3, cex.legend = 1.3, 
  edge.curved = 0.65, edge.label.cex = 1.3, edge.arrow.size = 0.8)

## ----code_ssplotMHMM, eval = FALSE, echo = TRUE---------------------------
#  mssplot(mhmm, ask = TRUE)

## ----setup, include=FALSE, cache=FALSE------------------------------------
library(knitr)
opts_chunk$set(concordance = TRUE, tidy = FALSE)
options(prompt = "R> ", continue = "+  ", width = 76, useFancyQuotes = FALSE)


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

## ----graphicalillustrations2, fig.width=6.5, fig.height=3.7, dev.args=list(pointsize=10), fig.keep='last', cache=FALSE, message=FALSE, echo=FALSE, fig.cap='Stacked sequence plot of the first ten individuals in the \\code{biofam} data plotted with the \\code{ssplot} function. The top plot shows the original sequences, and the three bottom plots show the sequences in the separate channels for the same individuals. The sequences are in the same order in each plot, i.e., the same row always matches the same individual.', fig.align='center'----
ssplot(list(biofam_seq[1:10,], marr_seq[1:10,], child_seq[1:10,],
  left_seq[1:10,]),
  sortv = "from.start", sort.channel = 1, type = "I",
  ylab = c("Original", "Marriage", "Parenthood", "Residence"),
  xtlab = 15:30, xlab = "Age", title = "Ten first sequences",
  title.n = FALSE, legend.prop = 0.63, ylab.pos = c(1, 1.5),
  ncol.legend = c(3, 1, 1, 1))

## ----plottingsequences, fig.width=5, fig.height=3, dev.args=list(pointsize=10), fig.cap="Stacked sequence plot of annual state distributions in the three-channel \\code{biofam} data. This is the default output of the \\code{ssplot} function. The labels for the channels are taken from the named list of state sequence objects, and the labels for the x axis ticks are taken from the column names of the first object.", fig.keep='last', fig.align='center', cache=FALSE, echo = FALSE----
ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq,
  "Residence" = left_seq))

## ----gridplot1, fig.width=5.5, fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Showing state distribution plots for women and men in the \\code{biofam} data. Two figures were defined with the \\code{ssp} function and then combined into one figure with the \\code{gridplot} function.", fig.align='center', fig.keep='last', cache = FALSE----
ssp_f <- ssp(list(marr_seq[biofam3c$covariates$sex == "woman",],
    child_seq[biofam3c$covariates$sex == "woman",],
    left_seq[biofam3c$covariates$sex == "woman",]),
  type = "I", sortv = "mds.obs", with.legend = FALSE, title = "Women", 
  ylab.pos = c(1, 2, 1), xtlab = 15:30, ylab = c("Married", "Children", 
    "Residence"))

ssp_m <- update(ssp_f, title = "Men", 
  x = list(marr_seq[biofam3c$covariates$sex == "man",],
    child_seq[biofam3c$covariates$sex == "man",],
    left_seq[biofam3c$covariates$sex == "man",]))

gridplot(list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
  legend.pos = "bottom", legend.pos2 = "top", row.prop = c(0.65, 0.35))

## ----code_mcHMM, cache=FALSE, echo = FALSE, message=FALSE, warning=TRUE, eval = TRUE----
mc_init <- c(0.9, 0.05, 0.02, 0.02, 0.01)

mc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0, 0.90, 0.05, 
  0.03, 0.02, 0, 0, 0.90, 0.07, 0.03, 0, 0, 0, 0.90, 0.10, 0, 0, 0,    
  0, 1), nrow = 5, ncol = 5, byrow = TRUE)

mc_emiss_marr <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 
  0.90, 0.05, 0.05, 0.90, 0.05, 0.30, 0.30, 0.40), nrow = 5, ncol = 3, 
  byrow = TRUE)

mc_emiss_child <- matrix(c(0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9,
  0.5, 0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_emiss_left <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9,
  0.5, 0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_obs <- list(marr_seq, child_seq, left_seq)

mc_emiss <- list(mc_emiss_marr, mc_emiss_child, mc_emiss_left)

mc_initmod <- build_hmm(observations = mc_obs, initial_probs = mc_init, 
  transition_probs = mc_trans, emission_probs = mc_emiss,
  channel_names = c("Marriage", "Parenthood", "Residence"))

# For CRAN vignette: load the estimated model object for speed-up
data("hmm_biofam")
# mc_fit <- fit_model(mc_initmod, em_step = FALSE, local_step = TRUE,
#   threads = 4)


## ----plottingHMM, out.width='\\linewidth', fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Illustrating a hidden Markov model as a directed graph. Pies represent five hidden states, with slices showing emission probabilities of combinations of observed states. States with emission probability less than 0.05 are combined into one slice. Edges show the transtion probabilities. Initial probabilities of hidden states are given below the pies.", fig.align='center', fig.keep='last', cache = FALSE----
plot(hmm_biofam, vertex.size = 50, vertex.label.dist = 1.5,
  edge.curved = c(0, 0.6, -0.8, 0.6, 0, 0.6, 0), legend.prop = 0.3, 
  combined.slice.label = "States with prob. < 0.05")

## ----graphicalillustrations5, out.width='\\linewidth', fig.height=3.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Another version of the hidden Markov model of Figure 4 with a different layout and modified labels, legends, and colors. All observed states are shown.", fig.align='center', fig.keep='last', cache = FALSE----
vertex_layout <- matrix(c(1, 2, 2, 3, 1, 0, 0.5, -0.5, 0, -1), 
  ncol = 2)

plot(hmm_biofam, layout = vertex_layout, xlim = c(0.5, 3.5), 
  ylim = c(-1.5, 1), rescale = FALSE, vertex.size = 50, 
  vertex.label.pos = c("left", "top", "bottom", "right", "left"),
  edge.curved = FALSE, edge.width = 1, edge.arrow.size = 1, 
  with.legend = "left", legend.prop = 0.4, label.signif = 1, 
  combine.slices = 0, cpal = colorpalette[[30]][c(14:5)])

## ----ssplotHMM, fig.width=5.5, fig.height=5.5, dev.args=list(pointsize=10), fig.cap="Using the \\code{ssplot} function for an \\code{hmm} object makes it easy to plot the observed sequences together with the most probable paths of hidden states given the model.", fig.align='center', fig.keep='last', cache = FALSE, echo = FALSE----
ssplot(hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden",
  title = "Observed and hidden state sequences", xtlab = 15:30, 
  xlab = "Age")

## ----code_settingdata, ref.label = 'settingdata', message=FALSE, warning = TRUE, echo = TRUE----
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

## ----code_plottingsequences, ref.label = 'plottingsequences', echo = TRUE, eval = FALSE----
#  ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq,
#    "Residence" = left_seq))

## ----code_graphicalillustrations2, ref.label = 'graphicalillustrations2', echo=TRUE, eval = FALSE----
#  ssplot(list(biofam_seq[1:10,], marr_seq[1:10,], child_seq[1:10,],
#    left_seq[1:10,]),
#    sortv = "from.start", sort.channel = 1, type = "I",
#    ylab = c("Original", "Marriage", "Parenthood", "Residence"),
#    xtlab = 15:30, xlab = "Age", title = "Ten first sequences",
#    title.n = FALSE, legend.prop = 0.63, ylab.pos = c(1, 1.5),
#    ncol.legend = c(3, 1, 1, 1))

## ----code_gridplot1, ref.label = 'gridplot1', echo=TRUE, eval = FALSE-----
#  ssp_f <- ssp(list(marr_seq[biofam3c$covariates$sex == "woman",],
#      child_seq[biofam3c$covariates$sex == "woman",],
#      left_seq[biofam3c$covariates$sex == "woman",]),
#    type = "I", sortv = "mds.obs", with.legend = FALSE, title = "Women",
#    ylab.pos = c(1, 2, 1), xtlab = 15:30, ylab = c("Married", "Children",
#      "Residence"))
#  
#  ssp_m <- update(ssp_f, title = "Men",
#    x = list(marr_seq[biofam3c$covariates$sex == "man",],
#      child_seq[biofam3c$covariates$sex == "man",],
#      left_seq[biofam3c$covariates$sex == "man",]))
#  
#  gridplot(list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
#    legend.pos = "bottom", legend.pos2 = "top", row.prop = c(0.65, 0.35))

## ----code_sc_buildHMM_random, cache=FALSE---------------------------------
sc_initmod_random <- build_hmm(observations = biofam_seq, n_states = 5)

## ----code_sc_initialvalues, cache=FALSE-----------------------------------
sc_init <- c(0.9, 0.06, 0.02, 0.01, 0.01)

sc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0.02, 0.80, 0.10, 
  0.05, 0.03, 0.02, 0.03, 0.80, 0.10, 0.05, 0.02, 0.03, 0.05, 0.80, 0.10, 
  0.02, 0.03, 0.05, 0.05, 0.85), nrow = 5, ncol = 5, byrow = TRUE)

sc_emiss <- matrix(NA, nrow = 5, ncol = 8)
sc_emiss[1,] <- seqstatf(biofam_seq[, 1:4])[, 2] + 0.1
sc_emiss[2,] <- seqstatf(biofam_seq[, 5:7])[, 2] + 0.1
sc_emiss[3,] <- seqstatf(biofam_seq[, 8:10])[, 2] + 0.1
sc_emiss[4,] <- seqstatf(biofam_seq[, 11:13])[, 2] + 0.1
sc_emiss[5,] <- seqstatf(biofam_seq[, 14:16])[, 2] + 0.1
sc_emiss <- sc_emiss / rowSums(sc_emiss)

rownames(sc_trans) <- colnames(sc_trans) <- rownames(sc_emiss) <-
  paste("State", 1:5)

colnames(sc_emiss) <- attr(biofam_seq, "labels")

sc_trans
round(sc_emiss, 3)

## ----code_sc_buildHMM, cache=FALSE----------------------------------------
sc_initmod <- build_hmm(observations = biofam_seq, initial_probs = sc_init,
  transition_probs = sc_trans, emission_probs = sc_emiss)

## ----code_sc_fitHMM, cache=FALSE------------------------------------------
sc_fit <- fit_model(sc_initmod)

## ----code_sc_results1, cache=FALSE----------------------------------------
sc_fit$logLik

## ----code_sc_results2, cache=FALSE----------------------------------------
sc_fit$model

## ----code_mcHMM2, echo = TRUE, message=TRUE, warnings=TRUE, eval = TRUE, cache = FALSE----
mc_init <- c(0.9, 0.05, 0.02, 0.02, 0.01)

mc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0, 0.90, 0.05, 0.03, 
  0.02, 0, 0, 0.90, 0.07, 0.03, 0, 0, 0, 0.90, 0.10, 0, 0, 0, 0, 1), 
  nrow = 5, ncol = 5, byrow = TRUE)

mc_emiss_marr <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 0.90, 
  0.05, 0.05, 0.90, 0.05, 0.30, 0.30, 0.40), nrow = 5, ncol = 3, 
  byrow = TRUE)

mc_emiss_child <- matrix(c(0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.5, 
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_emiss_left <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.5, 
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_obs <- list(marr_seq, child_seq, left_seq)

mc_emiss <- list(mc_emiss_marr, mc_emiss_child, mc_emiss_left)

mc_initmod <- build_hmm(observations = mc_obs, initial_probs = mc_init, 
  transition_probs = mc_trans, emission_probs = mc_emiss,
  channel_names = c("Marriage", "Parenthood", "Residence"))

# For CRAN vignette: load the estimated model object for speed-up
data("hmm_biofam")
# mc_fit <- fit_model(mc_initmod, em_step = FALSE, local_step = TRUE,
# threads = 4)


## ----code_mcHMM_BIC, cache=FALSE, echo = TRUE, message=FALSE, eval = TRUE----
# Vignette: already loaded hmm_biofam
# hmm_biofam <- mc_fit$model
BIC(hmm_biofam)

## ----code_MHMM, cache=FALSE, echo = TRUE, eval = TRUE, warning=TRUE-------
mc_init2 <- c(0.9, 0.05, 0.03, 0.02)

mc_trans2 <- matrix(c(0.85, 0.05, 0.05, 0.05, 0, 0.90, 0.05, 0.05, 0, 0, 
  0.95, 0.05, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

mc_emiss_marr2 <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 
  0.85, 0.10, 0.05, 0.80, 0.15), nrow = 4, ncol = 3, byrow = TRUE)

mc_emiss_child2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mc_emiss_left2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mhmm_init <- list(mc_init, mc_init2)

mhmm_trans <- list(mc_trans, mc_trans2)

mhmm_emiss <- list(list(mc_emiss_marr, mc_emiss_child, mc_emiss_left), 
  list(mc_emiss_marr2, mc_emiss_child2, mc_emiss_left2))

biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, 
  c(1908, 1935, 1945, 1957))
biofam3c$covariates$cohort <- factor(biofam3c$covariates$cohort, 
  labels=c("1909-1935", "1936-1945", "1946-1957"))

init_mhmm <- build_mhmm(observations = mc_obs, initial_probs = mhmm_init, 
  transition_probs = mhmm_trans, emission_probs = mhmm_emiss, 
  formula = ~sex + cohort, data = biofam3c$covariates, 
  channel_names = c("Marriage", "Parenthood", "Residence"), 
  cluster_names = c("Cluster 1", "Cluster 2"))

# vignette: less restarts and no parallelization
set.seed(1011)
mhmm_fit <- fit_model(init_mhmm, local_step = TRUE, threads = 1,
  control_em = list(restart = list(times = 10)))
mhmm <- mhmm_fit$model

## ----code_summaryMHMM, cache=FALSE----------------------------------------
summary(mhmm, conditional_se = FALSE)

## ----code_plottingHMMbasic, out.width='\\linewidth', fig.height=4, dev.args=list(pointsize=10), echo=TRUE, fig.align='center', fig.keep='last', cache = FALSE, eval = TRUE, fig.cap="A default plot of a hidden Markov model."----
plot(hmm_biofam)

## ----code_plottingHMM, ref.label='plottingHMM', echo=TRUE, eval = FALSE----
#  plot(hmm_biofam, vertex.size = 50, vertex.label.dist = 1.5,
#    edge.curved = c(0, 0.6, -0.8, 0.6, 0, 0.6, 0), legend.prop = 0.3,
#    combined.slice.label = "States with prob. < 0.05")

## ----code_graphicalillustrations5, ref.label = 'graphicalillustrations5', echo=TRUE, eval = FALSE----
#  vertex_layout <- matrix(c(1, 2, 2, 3, 1, 0, 0.5, -0.5, 0, -1),
#    ncol = 2)
#  
#  plot(hmm_biofam, layout = vertex_layout, xlim = c(0.5, 3.5),
#    ylim = c(-1.5, 1), rescale = FALSE, vertex.size = 50,
#    vertex.label.pos = c("left", "top", "bottom", "right", "left"),
#    edge.curved = FALSE, edge.width = 1, edge.arrow.size = 1,
#    with.legend = "left", legend.prop = 0.4, label.signif = 1,
#    combine.slices = 0, cpal = colorpalette[[30]][c(14:5)])

## ----code_ssplotHMM, ref.label = 'ssplotHMM', eval = FALSE, echo = TRUE----
#  ssplot(hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden",
#    title = "Observed and hidden state sequences", xtlab = 15:30,
#    xlab = "Age")

## ----code_plottingMHMMbasic, fig.width=6.5, fig.height=8, dev.args=list(pointsize=10), echo=TRUE, fig.align='center', fig.keep='last', cache = FALSE, eval = TRUE, fig.cap="Plotting submodels of an MHMM with the \\code{plot} method."----
plot(mhmm, interactive = FALSE, nrow = 2, legend.prop = 0.45,
  vertex.size = 50, vertex.label.cex = 1.3, cex.legend = 1.3, 
  edge.curved = 0.65, edge.label.cex = 1.3, edge.arrow.size = 0.8)

## ----code_ssplotMHMM, eval = FALSE, echo = TRUE---------------------------
#  mssplot(mhmm, ask = TRUE)

