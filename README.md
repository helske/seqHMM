[![Build Status](https://travis-ci.org/helske/seqHMM.svg?branch=master)](https://travis-ci.org/helske/seqHMM)

seqHMM: Hidden Markov Models for Life Sequences and Other Multivariate, Multichannel Categorical Time Series
====================================================================================================

Note! Function names and arguments were changed for version 0.12

The `seqHMM` package is designed for fitting hidden (or latent) Markov models (HMMs) and mixture hidden Markov models (MHMMs) for social sequence data and other categorical time series. The package supports models for one or multiple subjects with one or multiple interdependent sequences (channels). External covariates can be added to explain cluster membership in MHMMs. The package provides functions for evaluating and comparing models, as well as functions for easy plotting of multichannel sequences and hidden Markov models.

Maximum likelihood estimation via EM algorithm and direct numerical maximization with analytical gradients is supported. All main algorithms are written in C++.

Package is still under development and should be available at CRAN in 2015.

We would be happy to hear feedback! If you have any questions, comments, or wishes, please contact Satu Helske or Jouni Helske, firstname.lastname@jyu.fi.

If you want to try the `seqHMM` package, you can install it via the `devtools` package. Install also the latest version of `igraph`, which is needed for plotting HMMs. Just paste the following lines in R and you are ready to go.

```R
install.packages("devtools")
library(devtools)
install_github("helske/seqHMM")
install_github("igraph/rigraph")
```

Preview of the `seqHMM` package
---------------------------------------------------------------------------------

This example uses the `biofam` data from the `TraMineR` package. The data consist of a sample of 2000 individuals born between 1909 and 1972 constructed from the Swiss Household Panel (SHP) survey in 2002. The sequences consist of family life states from age 15 to 30 (in columns 10 to 25).


- 0 = "living with parents"
- 1 = "left home"
- 2 = "married"
- 3 = "left home+married"
- 4 = "child"
- 5 = "left home+child"
- 6 = "left home+married+child"
- 7 = "divorced"

For the functions of the `seqHMM` package, sequence data is given as a state sequence object (`stslist`) using the `seqdef` function in the `TraMineR` package. To show a more complex example, the original data is split into three separate channels. For the divorced state there is no information on children or residence, so these are assessed using the preceding states.

```
library(TraMineR)

data(biofam)

# Complete cases in sex, birthyear, and first nationality
bio <- biofam[complete.cases(biofam[c(2:4)]),]

# Sequence data for the first six individuals
head(bio[10:25])

# Building one channel per type of event (married, children, or left)
bf <- as.matrix(bio[, 10:25])
married <- bf == 2 | bf == 3 | bf == 6
children <-  bf == 4 | bf == 5 | bf == 6
left <- bf == 1 | bf == 3 | bf == 5 | bf == 6 | bf == 7

# Giving labels and modifying sequences

# Marriage
married[married == TRUE] <- "Married"
married[married == FALSE] <- "Single"
married[bf == 7] <- "Divorced"

# Parenthood
children[children == TRUE] <- "Children"
children[children == FALSE] <- "Childless"
# Divorced parents
div <- bf[
  (rowSums(bf == 7)>0 & rowSums(bf == 5)>0) | 
  (rowSums(bf == 7)>0 & rowSums(bf == 6)>0),
  ]
children[rownames(bf) %in% rownames(div) & bf == 7] <- "Children"

# Residence
left[left==TRUE] <- "Left home"
left[left==FALSE] <- "With parents"
# Divorced living with parents (before divorce)
wp <- bf[
  (rowSums(bf == 7) > 0 & rowSums(bf == 2) > 0 & rowSums(bf == 3) == 0 &
   rowSums(bf == 5) == 0 &  rowSums(bf == 6) == 0) |
  (rowSums(bf == 7) > 0 & rowSums(bf == 4) > 0 & rowSums(bf == 3) == 0 &
   rowSums(bf == 5) == 0 &  rowSums(bf==6) == 0),
  ]
left[rownames(bf) %in% rownames(wp) & bf == 7] <- "With parents"

# Building sequence objects (starting at age 15)
marr.seq <- seqdef(married, start = 15)
child.seq <- seqdef(children, start = 15)
left.seq <- seqdef(left, start = 15)

# Choosing colours for states
attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
```

### Plotting multichannel sequence data
Multichannel sequence data are easily plotted using the `ssplot` function (ssplot for Stacked Sequence Plot).

```
# Plotting state distribution plots of observations
ssplot(
  list(marr.seq, child.seq, left.seq), type = "d", plots = "obs", 
  title = "State distribution plots"
  )
```                  
![ssp1](https://github.com/helske/seqHMM/blob/master/Examples/ssp1.png)

Multiple `ssp` objects can also be plotted together in a grid.

```
# Preparing plots for women's state distributions
# Sorting by scores from multidimensional scaling
ssp_f2 <- ssp(
  list(marr.seq[bio$sex == "woman",], child.seq[bio$sex == "woman",],
       left.seq[bio$sex == "woman",]),
  type = "d", plots = "obs", border = NA, withlegend = FALSE,
  title = "State distributions for women", title.n = FALSE,
  ylab = c("Married", "Parenthood", "Left home"), ylab.pos = c(1,2,1),
  xlab = "Age", xtlab = 15:30
  )

# Same plot, but sequences instead of state distributions
ssp_f3 <- update(
  ssp_f2, type = "I", sortv = "mds.obs", title = "Sequences for women"
  )

# State distributions with men's data
ssp_m2 <- update(
  ssp_f2, x = list(marr.seq[bio$sex == "man",], child.seq[bio$sex == "man",], 
                   left.seq[bio$sex == "man",]),
  type = "d", plots = "obs", border = NA,
  title = "State distributions for men", title.n = FALSE,
  ylab = c("Married", "Parenthood", "Left home"), 
  withlegend = FALSE, ylab.pos = c(1,2,1)
  )

# Men's sequences
ssp_m3 <- update(
  ssp_m2, type = "I", sortv = "mds.obs", title = "Sequences for men"
  )

# Plotting state distributions and index plots of observations for women and men 
gridplot(list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), cols=2, byrow=TRUE, 
         row.prop=c(0.42,0.42,0.16))

```
![gridplot](https://github.com/helske/seqHMM/blob/master/Examples/gridplot.png)

### Fitting hidden Markov models

When fitting Hidden Markov models (HMMs), initial values for model parameters are first given to the `build_hmm` function. After that, the model is fitted with the `fit_hmm` function using EM algorithm, direct numerical estimation, or a combination of both.

```
# Initial values for emission matrices
B_marr <- matrix(NA, nrow=4, ncol=3)
B_marr[1,] <- seqstatf(marr.seq[, 1:4])[, 2] + 0.1
B_marr[2,] <- seqstatf(marr.seq[, 5:8])[, 2] + 0.1
B_marr[3,] <- seqstatf(marr.seq[, 9:12])[, 2] + 0.1
B_marr[4,] <- seqstatf(marr.seq[, 13:16])[, 2] + 0.1
B_marr <- B_marr / rowSums(B_marr)

B_child <- matrix(NA, nrow=4, ncol=2)
B_child[1,] <- seqstatf(child.seq[, 1:4])[, 2] + 0.1
B_child[2,] <- seqstatf(child.seq[, 5:8])[, 2] + 0.1
B_child[3,] <- seqstatf(child.seq[, 9:12])[, 2] + 0.1
B_child[4,] <- seqstatf(child.seq[, 13:16])[, 2] + 0.1
B_child <- B_child / rowSums(B_child)

B_left <- matrix(NA, nrow=4, ncol=2)
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
bHMM <- build_hmm(
  observations = list(marr.seq, child.seq, left.seq),
  initial_probs = initial_probs, transition_matrix = A, 
  emission_matrix = list(B_marr, B_child, B_left),
  channel_names = c("Marriage", "Parenthood", "Left home")
  )

# Fitting the HMM (using only the default MLSL algorithm)
HMM <- fit_hmm(bHMM)
HMM$logLik
# -14889.37

# Fitting with EM followed by MLSL algorithm
# Here leads to a better likelihood
HMM <- fit_hmm(bHMM, use_em = TRUE, use_nloptr = TRUE)
HMM$logLik
# -14883.86

```

### Plotting hidden Markov models

A simple `plot` method is used to show an `hmm` object as a graph. It shows hidden states as pie charts (vertices), with emission probabilities as slices and transition probabilities as arrows (edges). Initial probabilities are shown below the pies.

```
# Plot HMM
plot(HMM$model)
```
![HMMdefault](https://github.com/helske/seqHMM/blob/master/Examples/HMMdefault.png)

```

# A prettier version
plot(
  HMM$model,
  # larger vertices
  vertex.size = 45,
  # varying curvature of edges
  edge.curved = c(0,-0.7,0.6,0,-0.7,0),
  # legend with two columns and less space
  ncol.legend = 2, legend.prop = 0.4,
  # new label for combined slice
  combined.slice.label = "States with probability < 0.05")
```
![HMM](https://github.com/helske/seqHMM/blob/master/Examples/HMModel.png)

The `ssplot` function can also be used for plotting the observed states and/or the most probable paths of hidden states of a HMM.

```
# Plotting observations and hidden states
ssplot(HMM$model, plots = "both")
```
![sspboth_default](https://github.com/helske/seqHMM/blob/master/Examples/sspboth_default.png)
```
# Prettier version
ssplot(
  HMM$model, type="I", plots="both",
  # Sorting subjects according to multidimensional
  # scaling scores of the most probable hidden state paths
  sortv="mds.mpp", 
  # Naming the channels
  ylab=c("Children", "Married", "Left home"), 
  # Title for the plot
  title="Observed sequences and the 
most probable paths of hidden states",
  # Labels for hidden states (most common states)
  mpp.labels=c("1: Childless single, with parents", 
                "2: Childless single, left home",
                "3: Married without children",
                "4: Married parent, left home"),
  # Colours for hidden states
  mpp.col=c("olivedrab", "bisque", "plum", "indianred"),
  # Labels for x axis
  xtlab=15:30, xlab="Age",
  # Proportion for legends
  legend.prop=0.45)
```
![sspboth](https://github.com/helske/seqHMM/blob/master/Examples/sspboth.png)

### Computing likelihood and BIC

The `logLik` and `BIC` functions are used for model comparison with the log-likelihood or the Bayesian information criterion (BIC).

```
# Likelihood
logLik(HMM$model)
# -14883.86

# BIC
BIC(HMM$model)
# 30177.88
```

### Trimming HMMs

The `trim_hmm` function can be used to trim models by setting small probabilities to zero. Here the trimmed model led to model with slightly improved likelihood, so probabilities less than 0.01 could be set to zero.

```
trimmedHMM <- trim_hmm(HMM$model, maxit = 100, zerotol = 1e-04)
# "1 iteration(s) used."
# "Trimming improved log-likelihood, ll_trim-ll_orig = 5.57e-05"

# Emission probabilities of the original HMM
HMM$model$emiss
# $Marriage
#           symbol_names
# state_names     Divorced     Married       Single
#           1 0.000000e+00 2.67246e-20 1.000000e+00
#           2 2.306849e-45 8.27522e-12 1.000000e+00
#           3 4.978004e-02 9.50220e-01 1.046861e-12
#           4 1.802198e-02 9.49011e-01 3.296703e-02
# 
# $Parenthood
#           symbol_names
# state_names    Childless     Children
#           1 9.995067e-01 4.933051e-04
#           2 1.000000e+00 5.280791e-14
#           3 1.000000e+00 7.495200e-14
#           4 1.684766e-08 1.000000e+00
# 
# $`Left home`
#           symbol_names
# state_names    Left home With parents
#           1 2.064574e-21 1.000000e+00
#           2 1.000000e+00 1.730974e-15
#           3 7.024774e-01 2.975226e-01
#           4 1.000000e+00 8.557133e-53

# Emission probabilities of the trimmed HMM
trimmedHMM$emiss
# trimmedHMM$emiss
# $Marriage
#           symbol_names
# state_names   Divorced  Married     Single
#           1 0.00000000 0.000000 1.00000000
#           2 0.00000000 0.000000 1.00000000
#           3 0.04978004 0.950220 0.00000000
#           4 0.01802198 0.949011 0.03296703
# 
# $Parenthood
#           symbol_names
# state_names Childless     Children
#           1 0.9995067 0.0004933051
#           2 1.0000000 0.0000000000
#           3 1.0000000 0.0000000000
#           4 0.0000000 1.0000000000
# 
# $`Left home`
#           symbol_names
# state_names Left home With parents
#           1 0.0000000    1.0000000
#           2 1.0000000    0.0000000
#           3 0.7024774    0.2975226
#           4 1.0000000    0.0000000
```

### Converting multichannel to single channel models and data

The `mc_to_sc` function converts a multichannel model into a single channel representation. E.g. the `plot` function for `hmm` objects uses this type of conversion. The `seqHMM` package also includes a similar function `mc_to_sc_data` for merging multiple state sequence objects.

```
scHMM <- mc_to_sc(HMM$model)

ssplot(scHMM, plots = "both", sortv = "from.end", sort.channel = 0, 
       xtlab = 15:30, legend.prop = 0.5)
```
![scssp](https://github.com/helske/seqHMM/blob/master/Examples/scssp.png)

### Mixture hidden Markov models

A mixture hidden Markov model (MHMM) is, by definition, a mixture of HMMs that are fitted together. These are fitted and plotted with similar functions to ones presented before. Starting values are given as a list consisting of the parameter values for each cluster. The `buildMixHMM` function checks that the model is properly constructed before fitting with the `fitMixHMM`function. Trimming is called with the `trim_hmm`.
```
# Starting values for emission probabilities

# Cluster 1
alphabet(child.seq) # Checking for the order of observed states
B1_child <- matrix(c(0.99, 0.01, # High probability for childless
                     0.99, 0.01,
                     0.99, 0.01,
                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)

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
                    0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 2
B2_child <- matrix(c(0.99, 0.01, # High probability for childless
                     0.99, 0.01,
                     0.99, 0.01,
                     0.01, 0.99), nrow = 4, ncol = 2, byrow = TRUE)

B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
                    0.01, 0.01, 0.98,
                    0.01, 0.98, 0.01, # High probability for married
                    0.29, 0.7, 0.01), nrow = 4, ncol = 3, byrow = TRUE)

B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
                    0.99, 0.01,
                    0.99, 0.01,
                    0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 3
B3_child <- matrix(c(0.99, 0.01, # High probability for childless
                     0.99, 0.01,
                     0.01, 0.99,
                     0.99, 0.01,
                     0.01, 0.99,
                     0.01, 0.99), nrow = 6, ncol = 2, byrow = TRUE)

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
                    0.99, 0.01), nrow = 6, ncol = 2, byrow = TRUE)

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

# Creating covariate swiss
bio$swiss <- bio$nat_1_02 == "Switzerland"
bio$swiss[bio$swiss == TRUE] <- "Swiss"
bio$swiss[bio$swiss == FALSE] <- "Other"

# Build MHMM
bMHMM <- buildMixHMM(
  observations = list(marr.seq, child.seq, left.seq),
  transition_matrix = list(A1, A1, A2),
  emission_matrix = list(list(B1_marr, B1_child, B1_left), 
                        list(B2_marr, B2_child, B2_left),
                        list(B3_marr, B3_child, B3_left)),
  initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
  formula = ~ sex * birthyr + sex * swiss, data = bio, 
  cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
  channel_names = c("Marriage", "Parenthood", "Left home")
  )

MHMM <- fitMixHMM(bMHMM)

# Trim MHMM
trMHMM <- trim_hmm(MHMM$model, zerotol = 1e-04)
```

Parameter coefficients are stored in `beta`. The first cluster is the reference.
```
# trMHMM$beta
#                     Cluster 1    Cluster 2   Cluster 3
# (Intercept)                 0 -23.39241097 64.04751144
# sexwoman                    0  19.27426708 33.17879462
# birthyr                     0   0.01266753 -0.03279321
# swissSwiss                  0   0.11541402 -0.54893429
# sexwoman:birthyr            0  -0.01008529 -0.01725546
# sexwoman:swissSwiss         0   0.23815091  0.54652612
```
The most probable cluster for each individual is determined by the most probable path of hidden states. It is computed with the `mostProbablePath` function. The most probable clusters are stored in `mpp$cluster` and `mpp$classification_probabilities` gives the probability of each cluster (in columns) by the most probable cluster (rows). For individuals assigned to cluster 1, the average probability for cluster 1 is 0.84, 0.16 for cluster 2, and close to 0 for cluster 3. The highest probability for the assigned cluster is 0.94 for cluster 3.

```
# Computing most probable paths
mpp <- mostProbablePath(trMHMM)
# Assigning colours to hidden states
attr(mpp$mpp, "cpal") <- colorpalette[[14]]
# Number of individuals in each cluster
table(mpp$cluster)
# Cluster 1 Cluster 2 Cluster 3 
#       258      1236       281
# Cluster probabilities by the most probable cluster
mpp$classification_probabilities
#            Cluster 1  Cluster 2    Cluster 3
# Cluster 1 0.83738794 0.16213321 0.0004788443
# Cluster 2 0.07874202 0.87207334 0.0491846416
# Cluster 3 0.01403737 0.05066261 0.9353000274
```

### Plotting MHMMs

Also MHMMs are plotted with the `plot` function. The user can choose between an interactive mode (`interactive=TRUE`), where the model for each cluster is plotted separately, and a combined plot with all models at once.
```
# Plot mixture hidden Markov model
# Interactive plot, one cluster at a time
plot(trMHMM, interactive = TRUE)
```
![mixHMM1](https://github.com/helske/seqHMM/blob/master/Examples/mixHMM1.png)
![mixHMM2](https://github.com/helske/seqHMM/blob/master/Examples/mixHMM2.png)
![mixHMM3](https://github.com/helske/seqHMM/blob/master/Examples/mixHMM3.png)


```
# Plotting observed sequences and most probable hidden states
# Interactive plot, one cluster at a time
mssplot(
  trMHMM, plots = "both", sortv = "from.end", sort.channel = 1, 
  xtlab = 15:30, xlab = "Age"
  )

```
![mssplot1](https://github.com/helske/seqHMM/blob/master/Examples/mssplot1.png)
![mssplot2](https://github.com/helske/seqHMM/blob/master/Examples/mssplot2.png)
![mssplot3](https://github.com/helske/seqHMM/blob/master/Examples/mssplot3.png)




Coming later
---------------------------------------------------------------------------------------

<ul> 
 <li>Function for computing posterior probabilities</li>
 <li>Markov models</li>
 <li>Simulating sequences from HMMs</li>
</ul> 



