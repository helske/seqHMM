# Codes in the preview

library(TraMineR)

data(biofam3c)

# Building sequence objects (starting at age 15, with custom color palette)
marr.seq <- seqdef(
  biofam3c$married, start = 15, 
  cpal = c("#AB82FF", "#E6AB02", "#E7298A")
)
child.seq <- seqdef(
  biofam3c$children, start = 15,
  cpal = c("#66C2A5", "#FC8D62")
)
left.seq <- seqdef(
  biofam3c$left, start = 15,
  cpal = c("#A6CEE3", "#E31A1C")
)

# Plotting multichannel sequence data

# Plotting state distribution plots of observations
ssplot(
  list(marr.seq, child.seq, left.seq), type = "d", plots = "obs", 
  title = "State distribution plots")


# Preparing plots for women's state distributions
# Sorting by scores from multidimensional scaling
ssp_f2 <- ssp(
  list(marr.seq[bio$sex == "woman",], child.seq[bio$sex == "woman",],
       left.seq[bio$sex == "woman",]),
  type = "d", plots = "obs", border = NA, with.legend = FALSE,
  title = "State distributions for women", title.n = FALSE,
  ylab = c("Married", "Parenthood", "Left home"), ylab.pos = c(1,2,1),
  xlab = "Age", xtlab = 15:30)

# Same plot, but sequences instead of state distributions
ssp_f3 <- update(
  ssp_f2, type = "I", sortv = "mds.obs", title = "Sequences for women")

# State distributions with men's data
ssp_m2 <- update(
  ssp_f2, x = list(marr.seq[bio$sex == "man",], child.seq[bio$sex == "man",], 
                   left.seq[bio$sex == "man",]),
  title = "State distributions for men")

# Men's sequences
ssp_m3 <- update(
  ssp_m2, type = "I", sortv = "mds.obs", title = "Sequences for men")

# Plotting state distributions and index plots of observations for women and men 
gridplot(
  list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), cols = 2, byrow = TRUE, 
  row.prop = c(0.42, 0.42, 0.16))


# Fitting hidden Markov models

# Initial values for emission matrices
B_marr <- matrix(NA, nrow = 4, ncol = 3)
B_marr[1,] <- seqstatf(marr.seq[, 1:4])[, 2] + 0.1
B_marr[2,] <- seqstatf(marr.seq[, 5:8])[, 2] + 0.1
B_marr[3,] <- seqstatf(marr.seq[, 9:12])[, 2] + 0.1
B_marr[4,] <- seqstatf(marr.seq[, 13:16])[, 2] + 0.1
B_marr <- B_marr / rowSums(B_marr)

B_child <- matrix(NA, nrow = 4, ncol = 2)
B_child[1,] <- seqstatf(child.seq[, 1:4])[, 2] + 0.1
B_child[2,] <- seqstatf(child.seq[, 5:8])[, 2] + 0.1
B_child[3,] <- seqstatf(child.seq[, 9:12])[, 2] + 0.1
B_child[4,] <- seqstatf(child.seq[, 13:16])[, 2] + 0.1
B_child <- B_child / rowSums(B_child)

B_left <- matrix(NA, nrow = 4, ncol = 2)
B_left[1,] <- seqstatf(left.seq[, 1:4])[, 2] + 0.1
B_left[2,] <- seqstatf(left.seq[, 5:8])[, 2] + 0.1
B_left[3,] <- seqstatf(left.seq[, 9:12])[, 2] + 0.1
B_left[4,] <- seqstatf(left.seq[, 13:16])[, 2] + 0.1
B_left <- B_left / rowSums(B_left)

# Initial values for transition matrix
A <- matrix(c(0.9, 0.06, 0.03, 0.01,
              0,    0.9, 0.07, 0.03,
              0,      0,  0.9,  0.1,
              0,      0,    0,    1), nrow = 4, ncol = 4, byrow = TRUE)

# Initial values for initial state probabilities
initial_probs <- c(0.9, 0.07, 0.02, 0.01)

# Building the hidden Markov model with initial parameter values 
bhmm <- build_hmm(
  observations = list(marr.seq, child.seq, left.seq),
  initial_probs = initial_probs, transition_probs = A, 
  emission_probs = list(B_marr, B_child, B_left),
  channel_names = c("Marriage", "Parenthood", "Left home")
)

# Fitting the hmm (using only the default MLSL algorithm)
hmm <- fit_hmm(bhmm)
hmm$logLik

# Fitting with EM followed by MLSL algorithm
# Here leads to a better likelihood
hmm <- fit_hmm(bhmm, use_em = TRUE, use_nloptr = TRUE)
hmm$logLik


# Plot HMM
plot(hmm$model)

# A prettier version
plot(
  hmm$model,
  # larger vertices
  vertex.size = 45,
  # varying curvature of edges
  edge.curved = c(0,-0.7,0.6,0,-0.7,0),
  # legend with two columns and less space
  ncol.legend = 2, legend.prop = 0.4,
  # new label for combined slice
  combined.slice.label = "States with probability < 0.05")



# Plotting observations and hidden states
ssplot(hmm$model, plots = "both")

# Prettier version
ssplot(
  hmm$model, type = "I", plots = "both",
  # Sorting subjects according to multidimensional
  # scaling scores of the most probable hidden state paths
  sortv = "mds.mpp", 
  # Naming the channels
  ylab = c("Children", "Married", "Left home"), 
  # Title for the plot
  title = "Observed sequences and the 
  most probable paths of hidden states",
  # Labels for hidden states (most common states)
  mpp.labels = c("1: Childless single, with parents", 
                 "2: Childless single, left home",
                 "3: Married without children",
                 "4: Married parent, left home"),
  # Colours for hidden states
  mpp.col = c("olivedrab", "bisque", "plum", "indianred"),
  # Labels for x axis
  xtlab = 15:30, xlab = "Age",
  # Proportion for legends
  legend.prop = 0.45)


# Likelihood
logLik(hmm$model)

# BIC
BIC(hmm$model)


# Trimming HMM
trimmed_hmm <- trim_hmm(hmm$model, maxit = 100, zerotol = 1e-04)

# Emission probabilities of the original HMM
hmm$model$emiss

# Emission probabilities of the trimmed HMM
trimmed_hmm$emiss


# Multichannel to single channel HMM
sc_hmm <- mc_to_sc(hmm$model)

ssplot(
  sc_hmm, plots = "both", sortv = "from.end", sort.channel = 0, 
  xtlab = 15:30, legend.prop = 0.5)


# Mixture hidden Markov models
# Starting values for emission probabilities

# Cluster 1
alphabet(child.seq) # Checking for the order of observed states
B1_child <- matrix(c(0.99, 0.01, # High probability for childless
                     0.99, 0.01,
                     0.99, 0.01,
                     0.99, 0.01), 
                   nrow = 4, ncol = 2, byrow = TRUE)

alphabet(marr.seq)
B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
                    0.01, 0.01, 0.98,
                    0.01, 0.98, 0.01, # High probability for married
                    0.98, 0.01, 0.01), # High probability for divorced
                  nrow = 4, ncol = 3, byrow = TRUE)

alphabet(left.seq)
B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
                    0.99, 0.01, # High probability for having left home
                    0.99, 0.01,
                    0.99, 0.01), 
                  nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 2
B2_child <- matrix(c(0.99, 0.01, # High probability for childless
                     0.99, 0.01,
                     0.99, 0.01,
                     0.01, 0.99), 
                   nrow = 4, ncol = 2, byrow = TRUE)

B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
                    0.01, 0.01, 0.98,
                    0.01, 0.98, 0.01, # High probability for married
                    0.29, 0.7, 0.01), 
                  nrow = 4, ncol = 3, byrow = TRUE)

B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
                    0.99, 0.01,
                    0.99, 0.01,
                    0.99, 0.01), 
                  nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 3
B3_child <- matrix(c(0.99, 0.01, # High probability for childless
                     0.99, 0.01,
                     0.01, 0.99,
                     0.99, 0.01,
                     0.01, 0.99,
                     0.01, 0.99), 
                   nrow = 6, ncol = 2, byrow = TRUE)

B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
                    0.01, 0.01, 0.98,
                    0.01, 0.01, 0.98,
                    0.01, 0.98, 0.01, # High probability for married
                    0.01, 0.98, 0.01,
                    0.98, 0.01, 0.01), # High probability for divorced
                  nrow = 6, ncol = 3, byrow = TRUE)

B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
                    0.99, 0.01,
                    0.50, 0.50,
                    0.01, 0.99,
                    0.99, 0.01,
                    0.99, 0.01), 
                  nrow = 6, ncol = 2, byrow = TRUE)

# Starting values for transition matrices
A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
               0,    0.9, 0.07, 0.03,
               0,      0,  0.9,  0.1,
               0,      0,    0,    1),
             nrow = 4, ncol = 4, byrow = TRUE)

A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
               0,  0.7,  0.1,   0.1, 0.05, 0.05,
               0,    0, 0.85,  0.01,  0.1, 0.04,
               0,    0,    0,   0.9, 0.05, 0.05,
               0,    0,    0,     0,  0.9,  0.1,
               0,    0,    0,     0,    0,    1),
             nrow = 6, ncol = 6, byrow = TRUE)

# Starting values for initial state probabilities
initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)

# Build MHMM
bmhmm <- build_mhmm(
  observations = list(marr.seq, child.seq, left.seq),
  transition_probs = list(A1, A1, A2),
  emission_probs = list(list(B1_marr, B1_child, B1_left), 
                        list(B2_marr, B2_child, B2_left),
                        list(B3_marr, B3_child, B3_left)),
  initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
  formula = ~ sex * birthyr, data = biofam3c$covariates, 
  cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
  channel_names = c("Marriage", "Parenthood", "Left home")
)

mhmm <- fit_mhmm(bmhmm)

# Trim MHMM
tr_mhmm <- trim_hmm(mhmm$model, zerotol = 1e-04)

### Summary of MHMM

summ <- summary(tr_mhmm)
summ

### Plotting MHMMs

# Plot mixture hidden Markov model
# Interactive plot, one cluster at a time
plot(tr_mhmm, interactive = TRUE)

# Computing most probable paths
mpp <- hidden_paths(tr_mhmm)

# Plotting observed sequences and most probable hidden states
# Interactive plot, one cluster at a time
mssplot(
  tr_mhmm, plots = "both", sortv = "from.end", sort.channel = 1, 
  xtlab = 15:30, xlab = "Age"
)
