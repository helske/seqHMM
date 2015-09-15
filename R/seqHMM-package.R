#' The seqHMM package
#' 
#' The seqHMM package is designed for fitting hidden (or latent) Markov models (HMMs) and 
#' mixture hidden Markov models (MHMMs) for social sequence data and other categorical 
#' time series. The package supports models for one or multiple subjects with one or 
#' multiple interdependent sequences (channels). External covariates can be added to 
#' explain cluster membership in MHMMs. The package provides functions for evaluating 
#' and comparing models, as well as functions for easy plotting of multichannel sequences 
#' and hidden Markov models.
#' 
#' Maximum likelihood estimation via EM algorithm and direct numerical maximization 
#' with analytical gradients is supported. All main algorithms are written in C++.
#' 
#' @docType package
#' @name seqHMM
#' @aliases seqHMM
#' @examples 
#' require(TraMineR)
#' 
#' # Loading mvad and biofam3c data
#' data(mvad) 
#' data(biofam3c)
#' 
#' 
#' ###############################################################
#' 
#' ####### Dealing with multichannel sequence data #######
#' 
#' ##### Plotting multichannel data #####
#' 
#' # Creating sequence objects (data: biofam3c)
#' child.seq <- seqdef(biofam3c$children, start = 15)
#' marr.seq <- seqdef(biofam3c$married, start = 15)
#' left.seq <- seqdef(biofam3c$left, start = 15)
#' 
#' ## Choosing colors
#' attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#' 
#' # Plotting state distribution plots of observations
#' ssplot(list(child.seq, marr.seq, left.seq), type = "d", plots = "obs")
#' 
#' # Plotting sequence index plots of observations
#' ssplot(
#'   list(child.seq, marr.seq, left.seq), type = "I", plots = "obs",
#'   # Sorting subjects according to the beginning of the 2nd channel (marr.seq)
#'   sortv = "from.start", sort.channel = 2, 
#'   # Controlling the size, positions, and names for channel labels
#'   ylab.pos = c(1, 2, 1), cex.lab = 1, ylab = c("Children", "Married", "Left home"), 
#'   # Plotting without legend
#'   withlegend = FALSE
#'   )
#' 
#' ##### Plotting sequence data in a grid #####
#' 
#' # Creating sequence data (data: mvad)
#' 
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school", 
#' "training")
#' mvad.labels <- c("employment", "further education", "higher education", 
#'                  "joblessness", "school", "training")
#' mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, 
#'                    labels = mvad.labels, xtstep = 6)
#' 
#' # Defining plots for gridplot                   
#' ssp_m <- ssp(
#'   mvad.seq[mvad$male == "yes",], type = "d", withlegend = FALSE,
#'   title = "Men", ylab = NA, border = NA
#'   )
#' ssp_f <- update(ssp_m, x = mvad.seq[mvad$male == "no",], title = "Women")
#' 
#' # Plotting ssp_m and ssp_f in a grid
#' gridplot(list(ssp_m, ssp_f), ncol = 2, ncol.legend = 2)
#' 
#' 
#' ##### Converting multichannel to single-channel data #####
#' 
#' sc_biofam <- mc_to_sc_data(list(marr.seq, child.seq, left.seq))
#' 
#' ssplot(sc_biofam, type = "d", legend.prop = 0.5, ylab = "Proportion")
#' 
#' 
#' ###############################################################
#' 
#' ####### Building and fitting hidden Markov models (HMMs) #######
#' 
#' ##### Single-channel mvad data #####
#'
#' # Starting values for the emission matrix
#' emiss <- matrix(NA, nrow = 4, ncol = 6)
#' emiss[1,] <- seqstatf(mvad.seq[, 1:12])[, 2] + 0.1
#' emiss[2,] <- seqstatf(mvad.seq[, 13:24])[, 2] + 0.1
#' emiss[3,] <- seqstatf(mvad.seq[, 25:48])[, 2] + 0.1
#' emiss[4,] <- seqstatf(mvad.seq[, 49:70])[, 2] + 0.1
#' emiss <- emiss / rowSums(emiss)
#' 
#' # Starting values for the transition matrix
#' 
#' trans <-  matrix(c(0.80, 0.10, 0.05, 0.05,
#'                 0.05, 0.80, 0.10, 0.05,
#'                 0.05, 0.05, 0.80, 0.10,
#'                 0.05, 0.05, 0.10, 0.80), nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Starting values for initial state probabilities
#' initialpr <- c(0.3, 0.3, 0.2, 0.2)
#' 
#' # Building a hidden Markov model with starting values
#' bhmm_mvad <- build_hmm(
#'   observations = mvad.seq, transition_matrix = trans, 
#'   emission_matrix = emiss, initial_probs = initialpr
#' )
#' 
#' # Fitting with the EM algorithm
#' hmm_mvad <- fit_hmm(
#'   bhmm_mvad, em_step = TRUE, 
#'   global_step = FALSE, local_step = FALSE
#'   )
#'   
#'   
#' ##### Multichannel biofam3c data #####
#' 
#' # Starting values for emission matrices
#' B_marr <- matrix(NA, nrow=4, ncol=3)
#' B_marr[1,] <- seqstatf(marr.seq[, 1:4])[, 2] + 0.1
#' B_marr[2,] <- seqstatf(marr.seq[, 5:8])[, 2] + 0.1
#' B_marr[3,] <- seqstatf(marr.seq[, 9:12])[, 2] + 0.1
#' B_marr[4,] <- seqstatf(marr.seq[, 13:16])[, 2] + 0.1
#' B_marr <- B_marr / rowSums(B_marr)
#' 
#' B_child <- matrix(NA, nrow=4, ncol=2)
#' B_child[1,] <- seqstatf(child.seq[, 1:4])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 5:8])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 9:12])[, 2] + 0.1
#' B_child[4,] <- seqstatf(child.seq[, 13:16])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_left <- matrix(NA, nrow=4, ncol=2)
#' B_left[1,] <- seqstatf(left.seq[, 1:4])[, 2] + 0.1
#' B_left[2,] <- seqstatf(left.seq[, 5:8])[, 2] + 0.1
#' B_left[3,] <- seqstatf(left.seq[, 9:12])[, 2] + 0.1
#' B_left[4,] <- seqstatf(left.seq[, 13:16])[, 2] + 0.1
#' B_left <- B_left / rowSums(B_left)
#' 
#' # Initial values for transition matrix
#' A <- matrix(c(0.9, 0.06, 0.03, 0.01,
#'               0,    0.9, 0.07, 0.03,
#'               0,      0,  0.9,  0.1,
#'               0,      0,    0,    1), nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities
#' initial_probs <- c(0.9, 0.07, 0.02, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bhmm_biofam <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = initial_probs
#'   )
#' 
#' # Fitting with default steps:
#' # Step 1) EM algorithm
#' # Step 2) Global optimization via MLSL_LDS with LBFGS as local optimizer;
#' #         3000 evaluations, unlimited time
#' # Step 3) Local optimization with LBFGS algorithm for "final polishing";
#'           3000 evaluations, unlimited time
#' # Note: estimation time limited to 60 seconds by default
#' \dontrun{
#' hmm_biofam <- fit_hmm(
#'   bhmm_biofam, 
#'   control_global = list(maxeval = 3000, maxtime = 0),
#'   control_local = list(maxeval = 3000, maxtime = 0)
#'   )
#' }
#' 
#' 
#' ###############################################################
#' 
#' ####### Building and fitting mixture hidden Markov models (MHMMs) #######
#' 
#' ##### Single-channel mvad data #####
#' 
#' # Starting values for emission matrices
#' B1 <- matrix(c(0.26, 0.39, 0.01, 0.06, 0.04, 0.24,
#'                0.58, 0.12, 0.09, 0.10, 0.01, 0.10,
#'                0.73, 0.02, 0.09, 0.13, 0.01, 0.02), nrow = 3, ncol = 6, byrow = TRUE)
#' 
#' B2 <- matrix(c(0.01, 0.02, 0.01, 0.01, 0.94, 0.01,
#'                0.05, 0.06, 0.15, 0.01, 0.72, 0.01,
#'                0.19, 0.13, 0.60, 0.01, 0.05, 0.02,
#'                0.32, 0.03, 0.60, 0.03, 0.01, 0.01), nrow = 4, ncol = 6, byrow = TRUE)
#' 
#' # Starting values for transition matrices
#' 
#' A1 <-  matrix(c(0.80, 0.10, 0.10,
#'                 0.10, 0.80, 0.10,
#'                 0.10, 0.10, 0.80), nrow=3, ncol=3, byrow=TRUE)
#' 
#' A2 <-  matrix(c(0.80, 0.10, 0.05, 0.05,
#'                 0.05, 0.80, 0.10, 0.05,
#'                 0.05, 0.05, 0.80, 0.10,
#'                 0.05, 0.05, 0.10, 0.80), nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Starting values for initial state probabilities
#' initial_probs1 <- c(0.4, 0.3, 0.3)
#' initial_probs2 <- c(0.3, 0.3, 0.2, 0.2)
#' 
#' # Building a mixture hidden Markov model with the starting values
#' bmhmm_mvad <- build_mhmm(
#'   observations = mvad.seq, 
#'   transition_matrix = list(A1, A2), 
#'   emission_matrix = list(B1, B2), 
#'   initial_probs = list(initial_probs1, initial_probs2)
#' )
#' 
#' # Fitting a MHMM with the EM algorithm
#' mhmm_mvad <- fit_mhmm(
#'   bmhmm_mvad, em_step = TRUE, global_step = FALSE, local_step = FALSE
#'   )
#'   
#'   
#' ##### Multichannel biofam3c data #####
#' 
#' # Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' B1_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' alphabet(marr.seq)                      
#' B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                     nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01, # High probability for having left home
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' # Cluster 2
#' B2_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.01, 0.99), nrow = 4, ncol = 2, byrow = TRUE)
#'                      
#' B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.29, 0.7, 0.01),
#'                    nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE) 
#' 
#' # Cluster 3
#' B3_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.01, 0.99), nrow = 6, ncol = 2, byrow = TRUE)
#' 
#' B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                    nrow = 6, ncol = 3, byrow = TRUE)                   
#' 
#' B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.50, 0.50,
#'                     0.01, 0.99,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 6, ncol = 2, byrow = TRUE) 
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'                  0,    0.9, 0.07, 0.03, 
#'                  0,      0,  0.9,  0.1, 
#'                  0,      0,    0,    1), 
#'              nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
#'                  0,  0.7,  0.1,   0.1, 0.05, 0.05,
#'                  0,    0, 0.85,  0.01,  0.1, 0.04,
#'                  0,    0,    0,   0.9, 0.05, 0.05,
#'                  0,    0,    0,     0,  0.9,  0.1,
#'                  0,    0,    0,     0,    0,    1), 
#'              nrow = 6, ncol = 6, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Birth cohort
#' biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
#' biofam3c$covariates$cohort <- factor(
#'   biofam3c$covariates$cohort, labels=c("1909-1935", "1936-1945", "1946-1957")
#' )
#' 
#' # Build mixture HMM
#' bmhmm_biofam <- build_mhmm(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_matrix = list(A1,A1,A2),
#'   emission_matrix = list(list(B1_child, B1_marr, B1_left),
#'                         list(B2_child, B2_marr, B2_left), 
#'                         list(B3_child, B3_marr, B3_left)),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   formula = ~ sex + cohort, data = biofam3c$covariates,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Parenthood", "Marriage", "Left home")
#'   )
#' 
#' # Fitting the model with the EM algorithm
#' mhmm_biofam <- fit_mhmm(
#'   bmhmm_biofam, em_step = TRUE, global_step = FALSE, local_step = FALSE
#'   )
#' 
#' # Fitting with default steps
#' # Step 1) EM algorithm
#' # Step 2) Global optimization via MLSL_LDS with LBFGS as local optimizer;
#' #         3000 evaluations, unlimited time
#' # Step 3) Local optimization with LBFGS algorithm for "final polishing";
#'           3000 evaluations, unlimited time
#' # Note: estimation time limited to 60 seconds by default
#' \dontrun{
#' mhmm_biofam <- fit_mhmm(
#'   control_em = list(maxeval = 10), control_global = list(maxeval = 3000, maxtime = 0),
#'   control_local = list(maxeval = 3000, maxtime = 0)
#'   )
#' }
#' 
#' ###############################################################
#' 
#' ####### Model evaluation and manipulation #######
#' 
#' # Loading pre-fitted models (HMM and MHMM)
#' data(hmm_biofam)
#' data(mhmm_biofam)
#' 
#' 
#' ##### Log-likelihood and BIC for model comparison #####
#' logLik(hmm_biofam)
#' BIC(hmm_biofam)
#' 
#' logLik(mhmm_biofam)
#' BIC(mhmm_biofam)
#' 
#' 
#' ##### Trimming (mixture) hidden Markov models #####
#' # i.e. setting small (< zerotol) parameter values to zero
#' 
#' trhmm_biofam <- trim_hmm(hmm_biofam, zerotol = 1e-04)
#' trmhmm_biofam <- trim_hmm(mhmm_biofam, zerotol = 1e-03)
#' 
#' 
#' ##### Converting multichannel models to single-channel models #####
#' 
#' schmm_biofam <- mc_to_sc(hmm_biofam)
#' scmhmm_biofam <- mc_to_sc(mhmm_biofam) 
#' 
#' 
#' ##### Summary of the MHMM #####
#' 
#' summary(mhmm_biofam)
#' 
#' 
#' ###############################################################
#' 
#' ####### Plotting hidden Markov models #######
#' 
#' # Default plotting 
#' plot(hmm_biofam)
#' 
#' # Plotting HMM with
#' plot(hmm_biofam,
#'      # larger vertices
#'      vertex.size = 40,
#'      # varying curvature of edges
#'      edge.curved = c(0, 0.7, -0.6, 0, 0.7, 0),
#'      # legend with two columns and less space
#'      ncol.legend = 2, legend.prop = 0.4,
#'      # new label for combined slice
#'      combined.slice.label = "States with probability < 0.05",
#'      # new color palette (1 color for each combined observation in the data)
#'      cpal = colorpalette[[10]]
#'      )
#' 
#' 
#' ##### Plotting mixture hidden Markov models #####
#' 
#' # Plotting only the first cluster
#' plot(mhmm_biofam, which.plots = 1)
#' 
#' \dontrun{
#' # Plotting each cluster (change with Enter)
#' plot(mhmm_biofam)
#' 
#' # Choosing the cluster (one at a time)
#' # Exit with 0
#' plot(mhmm_biofam, ask = TRUE)
#' }
#' 
#' 
#' # Binomial regression
#'
#' require("MASS")
#' data("birthwt")
#' 
#' b1 <- matrix(c(1,0),ncol=2)
#' b2 <- matrix(c(0,1),ncol=2)
#' 
#' a <- matrix(1)
#' model <- build_mhmm(obs= seqdef(birthwt$low), low ~ age + lwt + smoke + ht, birthwt,
#'   transition_matrix = list(a, a), initial_probs = list(1, 1), emission_matrix = list(b1, b2))
#' fit <- fit_mhmm(model)
#' summary(fit$model)[c("beta", "beta_se", "logLik")]
#' summary(glm(low ~ age + lwt + smoke + ht, binomial, data = bwt))
#' 
# multinomial regression
#' 
#' require("nnet")
#' 
#' set.seed(123)
#' n <- 100
#' X <- cbind(1, x1 = runif(n, 0, 1), x2 =  runif(n, 0, 1))
#' beta <- cbind(0,c(-2, 5, -2), c(0, -2, 2))
#' pr <- exp(X %*% beta)  + rnorm(n*3)
#' pr <- pr/rowSums(pr)
#' y <- apply(pr, 1, which.max)
#' table(y)
#' 
#' a <- matrix(1)
#' b1 <- matrix(c(1, 0, 0),ncol=3)
#' b2 <- matrix(c(0, 1, 0),ncol=3)
#' b3 <- matrix(c(0, 0, 1),ncol=3)
#' model <- build_mhmm(obs= seqdef(y), ~ x1 + x2,  data.frame(X[, -1]),
#'   transition_matrix = list(a, a, a),initial_probs = list(1, 1, 1), emission_matrix = list(b1, b2, b3))
#' fit <- fit_mhmm(model, local_step = FALSE, global_step = FALSE)
#' summary(fit$model)[c("beta", "beta_se", "logLik")]
#' BIC(fit$model)
#' multinom(y ~ x1 + x2, data = data.frame(X[,-1]))
#' 
NULL



