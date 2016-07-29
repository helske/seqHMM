################## Testing gridplot ######################

library(seqHMM)

data(biofam3c)

# Creating sequence objects
child.seq <- seqdef(biofam3c$children, start = 15)
marr.seq <- seqdef(biofam3c$married, start = 15)
left.seq <- seqdef(biofam3c$left, start = 15)

## Choosing colors
attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")


# Preparing plot for state distribution plots of observations for women
ssp_f <- ssp(
  list(child.seq[biofam3c$covariates$sex == "woman",],
       marr.seq[biofam3c$covariates$sex == "woman",],
       left.seq[biofam3c$covariates$sex == "woman",]),
  type = "d", plots = "obs", title = "Women",
  ylab = c("Children", "Married", "Left home"))

# Preparing plot for state distribution plots of observations for men
ssp_m <- ssp(
  list(child.seq[biofam3c$covariates$sex == "man",],
       marr.seq[biofam3c$covariates$sex == "man",],
       left.seq[biofam3c$covariates$sex == "man",]),
  type = "d", plots = "obs", title = "Men",
  ylab = c("Children", "Married", "Left home"))

# Preparing plots for women's state distributions
ssp_f2 <- ssp(
  list(marr.seq[biofam3c$covariates$sex == "woman",],
       child.seq[biofam3c$covariates$sex == "woman",],
       left.seq[biofam3c$covariates$sex == "woman",]),
  type = "d", border = NA, withlegend = FALSE,
  title = "State distributions for women", title.n = FALSE, xtlab = 15:30,
  ylab.pos = c(1, 2, 1), ylab = c("Married", "Children", "Left home"))

# The same plot with sequences instead of state distributions
ssp_f3 <- update(
  ssp_f2, type = "I", sortv="mds.obs", title = "Sequences for women")

# State distributions with men's data
ssp_m2 <- update(
  ssp_f2, title = "State distributions for men",
  x = list(marr.seq[biofam3c$covariates$sex == "man",],
           child.seq[biofam3c$covariates$sex == "man",],
           left.seq[biofam3c$covariates$sex == "man",]))

# Men's sequences
ssp_m3 <- update(
  ssp_m2, type = "I", sortv="mds.obs", title = "Sequences for women")

data("hmm_biofam")
data("hmm_mvad")
data("mhmm_biofam")
sep_mhmm <- separate_mhmm(mhmm_biofam)

ssp_bf <- ssp(hmm_biofam, withlegend = FALSE)
ssp_mvad <- ssp(hmm_mvad)

mhmm_1 <- ssp(sep_mhmm[[1]], withlegend = FALSE)
mhmm_2 <- update(mhmm_1, x = sep_mhmm[[2]])
mhmm_3 <- update(mhmm_1, x = sep_mhmm[[3]])

############### Rows & columns ##################

gridplot(list(ssp_f2, ssp_m2))
gridplot(list(ssp_f2, ssp_m2), ncol = 2)
gridplot(list(ssp_f2, ssp_m2), nrow = 2)
gridplot(list(ssp_f2, ssp_m2), ncol = 2, nrow = 2)

gridplot(list(mhmm_1, mhmm_2, mhmm_3))
gridplot(list(mhmm_1, mhmm_2, mhmm_3), ncol = 2)
gridplot(list(mhmm_1, mhmm_2, mhmm_3), nrow = 3)
gridplot(list(mhmm_1, mhmm_2, mhmm_3), ncol = 2, nrow = 3)

############### Legend positions ################

# On right
gridplot(list(ssp_f2, ssp_m2), legend.pos = "right")
gridplot(list(ssp_f2, ssp_m2), ncol = 2, legend.pos = "right")
gridplot(list(ssp_f2, ssp_m2), nrow = 2, legend.pos = "right")
gridplot(list(ssp_f2, ssp_m2), ncol = 2, nrow = 2, legend.pos = "right")

gridplot(list(mhmm_1, mhmm_2, mhmm_3), legend.pos = "right",
         ncol.legend = 2, col.prop = c(0.4, 0.6))
gridplot(list(mhmm_1, mhmm_2, mhmm_3), ncol = 2, legend.pos = "right",
         ncol.legend = 2, col.prop = c(0.4, 0.6))
gridplot(list(mhmm_1, mhmm_2, mhmm_3), nrow = 3, legend.pos = "right",
         ncol.legend = 2, col.prop = c(0.4, 0.6))
gridplot(list(mhmm_1, mhmm_2, mhmm_3), ncol = 2, nrow = 3, legend.pos = "right",
         ncol.legend = 2, col.prop = c(0.4, 0.6))

# auto = bottom
gridplot(list(ssp_f2, ssp_m2), legend.pos = "auto")
gridplot(list(ssp_f2, ssp_m2), ncol = 2, legend.pos = "auto")
gridplot(list(ssp_f2, ssp_m2), nrow = 2, legend.pos = "auto")
gridplot(list(ssp_f2, ssp_m2), ncol = 2, nrow = 2, legend.pos = "auto")


################ Different data ###############

gridplot(list(ssp_bf, ssp_mvad)) # Expect warning
gridplot(list(ssp_bf, ssp_mvad), withlegend = FALSE)
