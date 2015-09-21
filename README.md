[![Build Status](https://travis-ci.org/helske/seqHMM.svg?branch=master)](https://travis-ci.org/helske/seqHMM)

### Note! Function names and arguments were changed in August and September 2015 (version 1.0.0 onwards). Sorry for the inconvenience.


seqHMM: Hidden Markov Models for Life Sequences and Other Multivariate, Multichannel Categorical Time Series
====================================================================================================

The `seqHMM` package is designed for fitting hidden (or latent) Markov models (HMMs) and mixture hidden Markov models (MHMMs) for social sequence data and other categorical time series. The package supports models for one or multiple subjects with one or multiple interdependent sequences (channels). External covariates can be added to explain cluster membership in MHMMs. The package provides functions for evaluating and comparing models, as well as functions for easy plotting of multichannel sequences and hidden Markov models.

Maximum likelihood estimation via EM algorithm and direct numerical maximization with analytical gradients is supported. All main algorithms are written in C++.

Package is still under development and should be available at CRAN in 2015.

We would be happy to hear feedback! If you have any questions, comments, or wishes, please contact Satu Helske or Jouni Helske, firstname.lastname@jyu.fi. You can also add a new issue here in GitHub.

If you want to try the `seqHMM` package, you can install it via the `devtools` package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/seqHMM")
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

For the functions of the `seqHMM` package, sequence data is given as a state sequence object (`stslist`) using the `seqdef` function in the `TraMineR` package. To show a more complex example, the original data is split into three separate channels. This data is pre-generated and stored as
`biofam3c`. It contains a list of three sequence data sets and a data frame including the covariates. Find more information and the code for the conversion by typing `help(biofam3c)`.

```
library(TraMineR)

data(biofam3c)

# Building sequence objects (starting at age 15)
marr.seq <- seqdef(biofam3c$married, start = 15)
child.seq <- seqdef(biofam3c$children, start = 15)
left.seq <- seqdef(biofam3c$left, start = 15)

# Choosing colours for states
attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
```

### Plotting multichannel sequence data
Multichannel sequence data are easily plotted using the `ssplot` function (ssplot for Stacked Sequence Plot). It can plot two types of plots: sequences themselves (index plots) or state distributions at each time point. Plotting index plots with large sequence data can take time, especially if the sequences are also sorted at the same time, so by default the function plots states distributions.

```
# Plotting state distribution plots of observations
ssplot(
  list("Marriage" = marr.seq, "Parenthood" = child.seq, "Residence" = left.seq), 
  title = "State distribution plots")
```                  
![ssp1](https://github.com/helske/seqHMM/blob/master/Examples/ssp1.png)

Multiple `ssp` objects can also be plotted together in a grid.

```
# Preparing plots for women's state distributions
# Sorting by scores from multidimensional scaling
ssp_f2 <- ssp(
  list(marr.seq[biofam3c$covariates$sex == "woman",], 
       child.seq[biofam3c$covariates$sex == "woman",],
       left.seq[biofam3c$covariates$sex == "woman",]),
  type = "d", plots = "obs", border = NA, withlegend = FALSE,
  title = "State distributions for women", title.n = FALSE,
  ylab = c("Married", "Parenthood", "Left home"), ylab.pos = c(1, 2, 1),
  xlab = "Age", xtlab = 15:30)

# Same plot, but sequences instead of state distributions
ssp_f3 <- update(
  ssp_f2, type = "I", sortv = "mds.obs", title = "Sequences for women")

# State distributions with men's data
ssp_m2 <- update(
  ssp_f2, x = list(marr.seq[biofam3c$covariates$sex == "man",], 
                   child.seq[biofam3c$covariates$sex == "man",], 
                   left.seq[biofam3c$covariates$sex == "man",]),
  title = "State distributions for men")

# Men's sequences
ssp_m3 <- update(
  ssp_m2, type = "I", sortv = "mds.obs", title = "Sequences for men")

# Plotting state distributions and index plots of observations for women and men 
gridplot(
  list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), ncol = 2, byrow = TRUE, 
  row.prop = c(0.42, 0.42, 0.16), legend.pos2 = "top")

```
![gridplot](https://github.com/helske/seqHMM/blob/master/Examples/gridplot.png)

### Fitting hidden Markov models

When fitting Hidden Markov models (HMMs), initial values for model parameters are first given to the `build_hmm` function. After that, the model is fitted with the `fit_hmm` function. The fitting function provides three estimation steps: 1) EM algorithm, 2) global optimization, and 3) local optimization. By default, all steps are performed. The results from a former step are used as starting values in a latter.

By default the `fit_hmm` function uses the multilevel single-linkage method (MLSL) with the LDS modification for global optimization. The MLSL method draws random starting points and performs a local optimization (LBFGS by default) from each. In order to reduce the computation time spent on non-global optima, the convergence tolerance of the local optimizer is set relatively large. At step 3, a local optimization (again LBFGS by default) is run with a lower tolerance to find the optimum with high precision.

There are some theoretical guarantees that the MLSL method used as the default optimizer in step 2 shoud ﬁnd all local optima in a ﬁnite number of local optimizations. Of course, it might not always succeed in a reasonable time. The EM algorithm can help in finding good boundaries for the search, especially with good starting values, but in some cases it can mislead. A good strategy is to try a couple of different fitting options with different combinations of the methods: e.g. all steps, only global and local steps, and a few evaluations of EM followed by global and local optimization.

By default, the estimation time is limited to 60 seconds in steps 2 and 3, so it is advisable to change the default settings for the final analysis. 

```
# Initial values for emission matrices
emiss_marr <- matrix(NA, nrow=4, ncol=3)
emiss_marr[1,] <- seqstatf(marr.seq[, 1:4])[, 2] + 0.1
emiss_marr[2,] <- seqstatf(marr.seq[, 5:8])[, 2] + 0.1
emiss_marr[3,] <- seqstatf(marr.seq[, 9:12])[, 2] + 0.1
emiss_marr[4,] <- seqstatf(marr.seq[, 13:16])[, 2] + 0.1
emiss_marr <- emiss_marr / rowSums(emiss_marr)

emiss_child <- matrix(NA, nrow=4, ncol=2)
emiss_child[1,] <- seqstatf(child.seq[, 1:4])[, 2] + 0.1
emiss_child[2,] <- seqstatf(child.seq[, 5:8])[, 2] + 0.1
emiss_child[3,] <- seqstatf(child.seq[, 9:12])[, 2] + 0.1
emiss_child[4,] <- seqstatf(child.seq[, 13:16])[, 2] + 0.1
emiss_child <- emiss_child / rowSums(emiss_child)

emiss_left <- matrix(NA, nrow=4, ncol=2)
emiss_left[1,] <- seqstatf(left.seq[, 1:4])[, 2] + 0.1
emiss_left[2,] <- seqstatf(left.seq[, 5:8])[, 2] + 0.1
emiss_left[3,] <- seqstatf(left.seq[, 9:12])[, 2] + 0.1
emiss_left[4,] <- seqstatf(left.seq[, 13:16])[, 2] + 0.1
emiss_left <- emiss_left / rowSums(emiss_left)

# Initial values for transition matrix
trans <- matrix(
  c(0.90, 0.06, 0.03, 0.01,
       0, 0.90, 0.07, 0.03,
       0,    0, 0.90, 0.10,
       0,    0,    0,    1), 
  nrow = 4, ncol = 4, byrow = TRUE)

# Initial values for initial state probabilities
initial_probs <- c(0.9, 0.07, 0.02, 0.01)

# Building the hidden Markov model with initial parameter values 
bhmm <- build_hmm(
  observations = list(marr.seq, child.seq, left.seq),
  initial_probs = initial_probs, transition_matrix = trans, 
  emission_matrix = list(emiss_marr, emiss_child, emiss_left),
  channel_names = c("Marriage", "Parenthood", "Residence"))

# Fitting the HMM in 3 estimation steps:
# step 1) EM algorithm
# step 2) global optimization (default: MLSL_LDS with LBFGS as local optimizer)
# step 3) local optimization (default: LBFGS) for "final polishing"
# Note: by default, estimation time limited to 60 seconds in steps 2 and 3
hmm <- fit_hmm(bhmm)
hmm$logLik
# -16854.16

# Only EM
hmm2 <- fit_hmm(bhmm, global_step = FALSE, local_step = FALSE)
hmm2$logLik
# -16854.16

# Only global optimization (3000 iterations, unlimited time)
hmm3 <- fit_hmm(bhmm, em_step = FALSE, local_step = FALSE,
control_global = list(maxeval = 3000, maxtime = 0))
hmm3$logLik
# -16856.78

```

### Plotting hidden Markov models

A simple `plot` method is used to show an `hmm` object as a graph. It shows hidden states as pie charts (vertices), with emission probabilities as slices and transition probabilities as arrows (edges). By default, initial probabilities are shown below the pies.

```
# Plot HMM
plot(hmm$model)
```
![HMMdefault](https://github.com/helske/seqHMM/blob/master/Examples/HMMdefault.png)

```

# A prettier version
plot(
  hmm$model,
  # larger vertices
  vertex.size = 45,
  # varying curvature of edges
  edge.curved = c(0, -0.7, 0.6, 0, -0.7, 0),
  # legend with two columns and less space
  ncol.legend = 2, legend.prop = 0.4,
  # new label for combined slice
  combined.slice.label = "States with probability < 0.05")
```
![HMM](https://github.com/helske/seqHMM/blob/master/Examples/HMModel.png)

The `ssplot` function can also be used for plotting the observed states and/or the most probable paths of hidden states of a HMM.

```
# Plotting observations and hidden states
ssplot(hmm$model, plots = "both", type = "I")
```
![sspboth_default](https://github.com/helske/seqHMM/blob/master/Examples/sspboth_default.png)
```
# Prettier version
ssplot(
  hmm$model, type = "I", plots = "both",
  # Sorting subjects according to multidimensional
  # scaling scores of the most probable hidden state paths
  sortv = "mds.hidden", 
  # Naming the channels
  ylab = c("Children", "Married", "Residence"), 
  # Title for the plot
  title = "Observed sequences and the 
most probable paths of hidden states",
  # Labels for hidden states (most common states)
  hidden.states.labels = c("1: Childless single, with parents", 
                           "2: Childless single, left home",
                           "3: Married without children",
                           "4: Married parent, left home"),
  # Colours for hidden states
  hidden.states.colors = c("olivedrab", "bisque", "plum", "indianred"),
  # Labels for x axis
  xtlab = 15:30, xlab = "Age",
  # Proportion for legends
  legend.prop = 0.45)
```
![sspboth](https://github.com/helske/seqHMM/blob/master/Examples/sspboth.png)

### Computing likelihood and BIC

The `logLik` and `BIC` functions are used for model comparison with the log-likelihood or the Bayesian information criterion (BIC).

```
# Likelihood
logLik(hmm$model)
# 'log Lik.' -16854.16 (df=25)

# BIC
BIC(hmm$model)
# 33967.66
```

### Trimming HMMs

The `trim_hmm` function can be used to trim models by setting small probabilities to zero. Here the trimmed model led to model with slightly improved likelihood, so probabilities less than 0.001 could be set to zero.

```
trimmed_hmm <- trim_hmm(hmm$model, maxit = 100, zerotol = 1e-03)
# "1 iteration(s) used."
# "Trimming improved log-likelihood, ll_trim-ll_orig = 4.28e-05"

# Transition probabilities of the original HMM
hmm$model$emission_matrix
# $Marriage
#            symbol_names
# state_names     divorced      married       single
#           1 0.000000e+00 3.026210e-22 1.000000e+00
#           2 1.793831e-50 2.301489e-09 1.000000e+00
#           3 4.713737e-02 9.528626e-01 6.419152e-14
#           4 1.747152e-02 9.497448e-01 3.278367e-02
# 
# $Parenthood
#            symbol_names
# state_names    childless     children
#           1 9.988180e-01 1.181965e-03
#           2 1.000000e+00 2.390547e-09
#           3 1.000000e+00 5.123136e-09
#           4 1.715309e-14 1.000000e+00
# 
# $`Residence`
#            symbol_names
# state_names    left home with parents
#           1 8.126227e-11 1.000000e+00
#           2 1.000000e+00 2.145277e-16
#           3 6.916852e-01 3.083148e-01
#           4 1.000000e+00 5.190911e-50

# Emission probabilities of the trimmed HMM
trimmed_hmm$emiss
# $Marriage
#            symbol_names
# state_names   divorced   married     single
#           1 0.00000000 0.0000000 1.00000000
#           2 0.00000000 0.0000000 1.00000000
#           3 0.04713737 0.9528626 0.00000000
#           4 0.01747154 0.9497448 0.03278367
# 
# $Parenthood
#            symbol_names
# state_names childless   children
#           1  0.998818 0.00118196
#           2  1.000000 0.00000000
#           3  1.000000 0.00000000
#           4  0.000000 1.00000000
# 
# $`Residence`
#            symbol_names
# state_names left home with parents
#           1 0.0000000    1.0000000
#           2 1.0000000    0.0000000
#           3 0.6916852    0.3083148
#           4 1.0000000    0.0000000
```

### Converting multichannel to single channel models and data

The `mc_to_sc` function converts a multichannel model into a single channel representation. E.g. the `plot` function for `hmm` objects uses this type of conversion. The `seqHMM` package also includes a similar function `mc_to_sc_data` for merging multiple state sequence objects.

```
sc_hmm <- mc_to_sc(hmm$model)

ssplot(sc_hmm, plots = "both", type = "I", sortv = "from.end", sort.channel = 0, 
       xtlab = 15:30, legend.prop = 0.5)
```
![scssp](https://github.com/helske/seqHMM/blob/master/Examples/scssp.png)

### Mixture hidden Markov models

A mixture hidden Markov model (MHMM) is, by definition, a mixture of HMMs that are fitted together. These are fitted and plotted with similar functions to ones presented before. Starting values are given as a list consisting of the parameter values for each cluster. The `build_mhmm` function checks that the model is properly constructed before fitting with the `fit_mhmm`function.
```
# Starting values for emission probabilities

# Cluster 1
alphabet(child.seq) # Checking for the order of observed states
emiss_1_child <- matrix(
  c(0.99, 0.01, # High probability for childless
    0.99, 0.01,
    0.99, 0.01,
    0.99, 0.01), 
  nrow = 4, ncol = 2, byrow = TRUE)

alphabet(marr.seq)
emiss_1_marr <- matrix(
  c(0.01, 0.01, 0.98, # High probability for single
    0.01, 0.01, 0.98,
    0.01, 0.98, 0.01, # High probability for married
    0.98, 0.01, 0.01), # High probability for divorced
  nrow = 4, ncol = 3, byrow = TRUE)

alphabet(left.seq)
emiss_1_left <- matrix(
  c(0.01, 0.99, # High probability for living with parents
    0.99, 0.01, # High probability for having left home
    0.99, 0.01,
    0.99, 0.01), 
  nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 2
emiss_2_child <- matrix(
  c(0.99, 0.01, # High probability for childless
    0.99, 0.01,
    0.99, 0.01,
    0.01, 0.99), 
  nrow = 4, ncol = 2, byrow = TRUE)

emiss_2_marr <- matrix(
  c(0.01, 0.01, 0.98, # High probability for single
    0.01, 0.01, 0.98,
    0.01, 0.98, 0.01, # High probability for married
    0.29, 0.7, 0.01), 
  nrow = 4, ncol = 3, byrow = TRUE)

emiss_2_left <- matrix(
  c(0.01, 0.99, # High probability for living with parents
    0.99, 0.01,
    0.99, 0.01,
    0.99, 0.01), 
  nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 3
emiss_3_child <- matrix(
  c(0.99, 0.01, # High probability for childless
    0.99, 0.01,
    0.01, 0.99,
    0.99, 0.01,
    0.01, 0.99,
    0.01, 0.99), 
  nrow = 6, ncol = 2, byrow = TRUE)

emiss_3_marr <- matrix(
  c(0.01, 0.01, 0.98, # High probability for single
    0.01, 0.01, 0.98,
    0.01, 0.01, 0.98,
    0.01, 0.98, 0.01, # High probability for married
    0.01, 0.98, 0.01,
    0.98, 0.01, 0.01), # High probability for divorced
  nrow = 6, ncol = 3, byrow = TRUE)

emiss_3_left <- matrix(
  c(0.01, 0.99, # High probability for living with parents
    0.99, 0.01,
    0.50, 0.50,
    0.01, 0.99,
    0.99, 0.01,
    0.99, 0.01), 
  nrow = 6, ncol = 2, byrow = TRUE)

# Starting values for transition matrices
trans_1 <- matrix(
  c(0.80, 0.16, 0.03, 0.01,
       0, 0.90, 0.07, 0.03,
       0,    0, 0.90, 0.10,
       0,    0,    0,    1),
  nrow = 4, ncol = 4, byrow = TRUE)

trans_2 <- matrix(
  c(0.80, 0.10, 0.05, 0.03, 0.01, 0.01,
       0, 0.70, 0.10, 0.10, 0.05, 0.05,
       0,    0, 0.85, 0.01, 0.10, 0.04,
       0,    0,    0, 0.90, 0.05, 0.05,
       0,    0,    0,    0, 0.90,  0.1,
       0,    0,    0,    0,    0,    1),
  nrow = 6, ncol = 6, byrow = TRUE)

# Starting values for initial state probabilities
initial_probs_1 <- c(0.9, 0.07, 0.02, 0.01)
initial_probs_2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)

# Birth cohort
biofam3c$covariates$cohort <- cut(
  biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
biofam3c$covariates$cohort <- factor(
  biofam3c$covariates$cohort, labels=c("1909-1935", "1936-1945", "1946-1957")
  )

# Build MHMM
bmhmm <- build_mhmm(
  observations = list(marr.seq, child.seq, left.seq),
  transition_matrix = list(trans_1, trans_1, trans_2),
  emission_matrix = list(list(emiss_1_marr, emiss_1_child, emiss_1_left), 
                        list(emiss_2_marr, emiss_2_child, emiss_2_left),
                        list(emiss_3_marr, emiss_3_child, emiss_3_left)),
  initial_probs = list(initial_probs_1, initial_probs_1, initial_probs_2),
  formula = ~sex * cohort, data = biofam3c$covariates, 
  cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
  channel_names = c("Marriage", "Parenthood", "Left home")
  )

mhmm <- fit_mhmm(bmhmm)
```

### Summary of MHMM

The `summary` method computes summaries of the MHMM, e.g. standard errors for covariates and prior and posterior probabilities for subjects. A `print` method shows some summaries of these: estimates and standard errors for covariates, log-likelihood and BIC, information on most probable clusters and prior and posterior probabilities. Parameter estimates for transitions, emissions, and initial probabilities are omitted by default. 

The classification table shows the mean probabilities of belonging to each cluster by the most probable cluster. The most probable cluster is determined by the posterior probabilities (or the most probable hidden state paths given by the `hidden_paths` function). A good model shoud have high proportions in the diagonal. Here, for individuals assigned to cluster 1, the average probability for cluster 1 is 0.84, 0.16 for cluster 2, and close to 0 for cluster 3. The highest probability for the assigned cluster is 0.93 for cluster 3.

```
summ_mhmm <- summary(mhmm$model)
names(summ_mhmm)
#  [1] "logLik"                          "BIC"                            
#  [3] "most_probable_cluster"           "coefficients"                   
#  [5] "coef_se"                         "prior_cluster_probabilities"    
#  [7] "posterior_cluster_probabilities" "classification_table"           
#  [9] "digits"                          "model" 
summ_mhmm
# Covariate effects :
# 
# Cluster 1 is the reference.
# 
# Cluster 2 :
#                           Estimate  Std. error
# (Intercept)                 1.1400       0.176
# sexwoman                   -0.2150       0.241
# cohort1936-1945             0.0829       0.239
# cohort1946-1957             0.1420       0.218
# sexwoman:cohort1936-1945    0.2960       0.329
# sexwoman:cohort1946-1957    0.0715       0.295
# 
# Cluster 3 :
#                           Estimate  Std. error
# (Intercept)                  0.430       0.197
# sexwoman                     0.149       0.263
# cohort1936-1945             -0.647       0.290
# cohort1946-1957             -0.899       0.269
# sexwoman:cohort1936-1945     0.212       0.387
# sexwoman:cohort1946-1957    -0.122       0.356
# 
# 
# Log-likelihood :
# 
# -12712.65
# 
# 
# BIC :
# 
# 26711.6
# 
# Means of prior cluster probabilities :
# 
# Cluster 1 Cluster 2 Cluster 3 
#     0.191     0.622     0.187
#
# 
# Most probable clusters :
# 
#             Cluster 1  Cluster 2  Cluster 3
# count             310       1364        326
# proportion      0.155      0.682      0.163
# 
# 
# Classification table :
# Mean cluster probabilities (in columns) by the most probable cluster (rows)
# 
#           Cluster 1 Cluster 2 Cluster 3
# Cluster 1     0.839     0.161     0.000
# Cluster 2     0.085     0.863     0.052
# Cluster 3     0.017     0.050     0.933
```

### Plotting MHMMs

Also MHMMs are plotted with the `plot` function. The user can choose between an interactive mode (`interactive=TRUE`), where the model for each cluster is plotted separately, and a combined plot with all models at once.
```
# Plot mixture hidden Markov model
# Interactive plot, one cluster at a time
plot(mhmm$model, interactive = TRUE)
```
![mixHMM1](https://github.com/helske/seqHMM/blob/master/Examples/mixHMM1.png)
![mixHMM2](https://github.com/helske/seqHMM/blob/master/Examples/mixHMM2.png)
![mixHMM3](https://github.com/helske/seqHMM/blob/master/Examples/mixHMM3.png)


```
# Plotting observed sequences and most probable hidden states
# Interactive plot, one cluster at a time
mssplot(
  mhmm$model, plots = "both", type = "I", sortv = "from.end", sort.channel = 1, 
  xtlab = 15:30, xlab = "Age")

```
![mssplot1](https://github.com/helske/seqHMM/blob/master/Examples/mssplot1.png)
![mssplot2](https://github.com/helske/seqHMM/blob/master/Examples/mssplot2.png)
![mssplot3](https://github.com/helske/seqHMM/blob/master/Examples/mssplot3.png)




Coming later
---------------------------------------------------------------------------------------

<ul>
 <li>Markov models</li>
</ul> 



