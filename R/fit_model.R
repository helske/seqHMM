#' Estimate Parameters of (Mixture) Hidden
#' Markov Models and Their Restricted Variants
#'
#' Function \code{fit_model} estimates the parameters of mixture hidden
#' Markov models and its restricted variants using maximimum likelihood.
#' Initial values for estimation are taken from the corresponding components
#' of the model with preservation of original zero probabilities.
#'
#' @export
#' @param model An object of class \code{hmm} or \code{mhmm}.
#' @param em_step Logical. Whether or not to use the EM algorithm at the start 
#'   of the parameter estimation. The default is \code{TRUE}.
#' @param global_step Logical. Whether or not to use global optimization via
#'   \code{\link{nloptr}} (possibly after the EM step). The default is \code{FALSE}.
#'@param local_step Logical. Whether or not to use local optimization via
#'   \code{\link{nloptr}} (possibly after the EM and/or global steps). The default is \code{FALSE}.
#' @param control_em Optional list of control parameters for the EM algorithm.
#'   Possible arguments are \describe{
#'   \item{maxeval}{The maximum number of iterations, the default is 1000.}
#'   \item{print_level}{The level of printing. Possible values are 0
#'   (prints nothing), 1 (prints information at the start and the end of the algorithm),
#'   2 (prints at every iteration),
#'   and for mixture models 3 (print also during optimization of coefficients).}
#'   \item{reltol}{Relative tolerance for convergence defined as
#'   \eqn{(logLik_new - logLik_old)/(abs(logLik_old) + 0.1)}.
#'   The default is 1e-12.}
#'   \item{restart}{A list containing options for possible EM restarts with the 
#'     following components:
#'   \describe{
#'   \item{times}{Number of restarts of the EM algorithm using random initial values. The default is 0, i.e. no restarts. }
#'   \item{transition}{Logical. Should the original transition probabilities be varied? The default is \code{TRUE}. }
#'   \item{emission}{Logical. Should the original emission probabilities be varied? The default is \code{TRUE}. }
#'   \item{sd}{Standard deviation for \code{rnorm} used in randomization. The default is 0.25.}
#'   \item{maxeval}{Maximum number of iterations, the default is 100.}
#'   \item{print_level}{Level of printing in restarted EM steps. The default is \code{control_em$print_level}. }
#'   \item{reltol}{Relative tolerance for convergence at restarted EM steps. The default is 1e-8.
#'   If the relative change of the final model of the restart phase is larger than the tolerance
#'   for the original EM phase, the final model is re-estimated with the original \code{reltol}
#'   and \code{maxeval} at the end of the EM step.}
#'   \item{n_optimum}{Save the log-likelihood values of the \code{n_optimum} best
#'   models (from all estimated models including the original).
#'   The default is \code{min(times + 1, 25)}.}
#'   }
#'   }
#'   }
#' @param control_global Optional list of additional arguments for
#'   \code{\link{nloptr}} argument \code{opts}. The default values are
#'   \describe{
#'    \item{algorithm}{\code{"NLOPT_GD_MLSL_LDS"}}
#'    \item{local_opts}{\code{list(algorithm = "NLOPT_LD_LBFGS", ftol_rel = 1e-6, xtol_rel = 1e-4)}}
#'    \item{maxeval}{\code{10000} (maximum number of iterations in global optimization algorithm.)}
#'    \item{maxtime}{\code{60} (maximum time for global optimization. Set to 0 for unlimited time.)}
#'}
#' @param lb,ub Lower and upper bounds for parameters in Softmax parameterization.
#' The default interval is \eqn{[pmin(-25, 2*initialvalues), pmax(25, 2*initialvalues)]}, 
#' except for beta coefficients,
#' where the scale of covariates is taken into account.
#' Note that it might still be a good idea to scale covariates around unit scale.
#' Bounds are used only in the global optimization step.
#'
#' @param control_local Optional list of additional arguments for
#'   \code{\link{nloptr}} argument \code{opts}. The default values are
#'   \describe{
#'    \item{algorithm}{\code{"NLOPT_LD_LBFGS"}}
#'    \item{ftol_rel}{\code{1e-10}}
#'    \item{xtol_rel}{\code{1e-8}}
#'    \item{maxeval}{\code{10000} (maximum number of iterations)}
#'   }
#' @param threads Number of threads to use in parallel computing. The default is 1.
#' @param log_space Make computations using log-space instead of scaling for greater
#' numerical stability at a cost of decreased computational performance. The default is \code{FALSE}.
#' @param ... Additional arguments to \code{nloptr}.
#' @return \describe{
#'   \item{logLik}{Log-likelihood of the estimated model. }
#'   \item{em_results}{Results after the EM step: log-likelihood (\code{logLik}), number of iterations
#'   (\code{iterations}), relative change in log-likelihoods between the last two iterations (\code{change}), and
#'   the log-likelihoods of the \code{n_optimum} best models after the EM step (\code{best_opt_restart}). }
#'   \item{global_results}{Results after the global step. }
#'   \item{local_results}{Results after the local step. }
#'   \item{call}{The matched function call. }
#'   }
#' @seealso \code{\link{build_hmm}},  \code{\link{build_mhmm}},
#' \code{\link{build_mm}},  \code{\link{build_mmm}}, and  \code{\link{build_lcm}}
#'   for constructing different types of models; \code{\link{summary.mhmm}}
#'   for a summary of a MHMM; \code{\link{separate_mhmm}} for reorganizing a MHMM into
#'   a list of separate hidden Markov models; \code{\link{plot.hmm}} and \code{\link{plot.mhmm}}
#'   for plotting model objects; and \code{\link{ssplot}} and \code{\link{mssplot}} for plotting
#'   stacked sequence plots of \code{hmm} and \code{mhmm} objects.
#' @details The fitting function provides three estimation steps: 1) EM algorithm,
#'   2) global optimization, and 3) local optimization. The user can call for one method
#'   or any combination of these steps, but should note that they are preformed in the
#'   above-mentioned order. The results from a former step are used as starting values
#'   in a latter, except for some of global optimization algorithms (such as MLSL and StoGO)
#'   which only use initial values for setting up the boundaries for the optimization.
#'
#'   It is possible to rerun the EM algorithm automatically using random starting
#'   values based on the first run of EM. Number of restarts is defined by
#'   the \code{restart} argument in \code{control_em}. As the EM algorithm is 
#'   relatively fast, this method might be preferred option compared to the proper 
#'   global optimization strategy of step 2.
#'
#'   The default global optimization method (triggered via \code{global_step = TRUE}) is
#'   the multilevel single-linkage method (MLSL) with the LDS modification (\code{NLOPT_GD_MLSL_LDS} as
#'   \code{algorithm} in \code{control_global}), with L-BFGS as the local optimizer.
#'   The MLSL method draws random starting points and performs a local optimization
#'   from each. The LDS modification uses low-discrepancy sequences instead of
#'   pseudo-random numbers as starting points and should improve the convergence rate.
#'   In order to reduce the computation time spent on non-global optima, the
#'   convergence tolerance of the local optimizer is set relatively large. At step 3,
#'   a local optimization (L-BFGS by default) is run with a lower tolerance to find the
#'   optimum with high precision.
#'
#'   There are some theoretical guarantees that the MLSL method used as the default
#'   optimizer in step 2 shoud find all local optima in a finite number of local
#'   optimizations. Of course, it might not always succeed in a reasonable time.
#'   The EM algorithm can help in finding good boundaries for the search, especially
#'   with good starting values, but in some cases it can mislead. A good strategy is to
#'   try a couple of different fitting options with different combinations of the methods:
#'   e.g. all steps, only global and local steps, and a few evaluations of EM followed by
#'   global and local optimization.
#'
#'   By default, the estimation time is limited to 60 seconds in global optimization step, so it is
#'   advisable to change the default settings for the proper global optimization.
#'
#'   Any algorithm available in the \code{nloptr} function can be used for the global and
#'   local steps.
#'
#' @examples
#'
#' # Hidden Markov model
#'
#' data("mvad", package = "TraMineR")
#'
#' mvad_alphabet <-
#'   c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad_labels <- c("employment", "further education", "higher education",
#'   "joblessness", "school", "training")
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet,
#'   states = mvad_scodes, labels = mvad_labels, xtstep = 6)
#'
#' attr(mvad_seq, "cpal") <- colorpalette[[6]]
#'
#' # Starting values for the emission matrix
#' emiss <- matrix(
#'   c(0.05, 0.05, 0.05, 0.05, 0.75, 0.05, # SC
#'     0.05, 0.75, 0.05, 0.05, 0.05, 0.05, # FE
#'     0.05, 0.05, 0.05, 0.4,  0.05, 0.4,  # JL, TR
#'     0.05, 0.05, 0.75, 0.05, 0.05, 0.05, # HE
#'     0.75, 0.05, 0.05, 0.05, 0.05, 0.05),# EM
#'   nrow = 5, ncol = 6, byrow = TRUE)
#'
#' # Starting values for the transition matrix
#' trans <- matrix(0.025, 5, 5)
#' diag(trans) <- 0.9
#'
#' # Starting values for initial state probabilities
#' initial_probs <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#'
#' # Building a hidden Markov model
#' init_hmm_mvad <- build_hmm(observations = mvad_seq,
#'   transition_probs = trans, emission_probs = emiss,
#'   initial_probs = initial_probs)
#'
#' \dontrun{
#' set.seed(21)
#' fit_hmm_mvad <- fit_model(init_hmm_mvad, control_em = list(restart = list(times = 100)))
#' hmm_mvad <- fit_hmm_mvad$model
#' }
#'
#' # save time, load the previously estimated model
#' data("hmm_mvad")
#'
#' # Markov model
#' 
#' set.seed(123)
#' init_mm_mvad <- build_mm(observations = mvad_seq,
#'   transition_probs = simulate_transition_probs(6),
#'   initial_probs = rep(1/6, 6))
#'
#' mm_mvad <- fit_model(init_mm_mvad)$model
#'
#' # Comparing likelihoods, MM fits better
#' logLik(hmm_mvad)
#' logLik(mm_mvad)
#'
#' \dontrun{
#' require("igraph") #for layout_in_circle
#'
#' plot(mm_mvad, layout = layout_in_circle, legend.prop = 0.3,
#'   edge.curved = 0.3, edge.label = NA,
#'   vertex.label.pos = c(0, 0, pi, pi, pi, 0))
#'
#' ##############################################################
#'
#'
#' #' # Three-state three-channel hidden Markov model
#' # See ?hmm_biofam for five-state version
#'
#' data("biofam3c")
#'
#' # Building sequence objects
#' marr_seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child_seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left_seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' # Define colors
#' attr(marr_seq, "cpal") <- c("violetred2", "darkgoldenrod2", "darkmagenta")
#' attr(child_seq, "cpal") <- c("darkseagreen1", "coral3")
#' attr(left_seq, "cpal") <- c("lightblue", "red3")
#'
#' # Starting values for emission matrices
#'
#' emiss_marr <- matrix(NA, nrow = 3, ncol = 3)
#' emiss_marr[1,] <- seqstatf(marr_seq[, 1:5])[, 2] + 1
#' emiss_marr[2,] <- seqstatf(marr_seq[, 6:10])[, 2] + 1
#' emiss_marr[3,] <- seqstatf(marr_seq[, 11:16])[, 2] + 1
#' emiss_marr <- emiss_marr / rowSums(emiss_marr)
#'
#' emiss_child <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_child[1,] <- seqstatf(child_seq[, 1:5])[, 2] + 1
#' emiss_child[2,] <- seqstatf(child_seq[, 6:10])[, 2] + 1
#' emiss_child[3,] <- seqstatf(child_seq[, 11:16])[, 2] + 1
#' emiss_child <- emiss_child / rowSums(emiss_child)
#'
#' emiss_left <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_left[1,] <- seqstatf(left_seq[, 1:5])[, 2] + 1
#' emiss_left[2,] <- seqstatf(left_seq[, 6:10])[, 2] + 1
#' emiss_left[3,] <- seqstatf(left_seq[, 11:16])[, 2] + 1
#' emiss_left <- emiss_left / rowSums(emiss_left)
#'
#' # Starting values for transition matrix
#' trans <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#'
#' # Starting values for initial state probabilities
#' inits <- c(0.9, 0.09, 0.01)
#'
#' # Building hidden Markov model with initial parameter values
#' init_hmm_bf <- build_hmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   transition_probs = trans,
#'   emission_probs = list(emiss_marr, emiss_child, emiss_left),
#'   initial_probs = inits)
#'
#' # Fitting the model with different optimization schemes
#'
#' # Only EM with default values
#' hmm_1 <- fit_model(init_hmm_bf)
#' hmm_1$logLik # -24179.1
#'
#' # Only L-BFGS
#' hmm_2 <- fit_model(init_hmm_bf, em_step = FALSE, local_step = TRUE)
#' hmm_2$logLik # -22267.75
#'
#' # Global optimization via MLSL_LDS with L-BFGS as local optimizer and final polisher
#' # This can be slow, use parallel computing by adjusting threads argument
#' # (here threads = 1 for portability issues)
#' hmm_3 <- fit_model(
#'   init_hmm_bf, em_step = FALSE, global_step = TRUE, local_step = TRUE,
#'   control_global = list(maxeval = 5000, maxtime = 0), threads = 1)
#' hmm_3$logLik # -21675.42
#'
#' # EM with restarts, much faster than MLSL
#' set.seed(123)
#' hmm_4 <- fit_model(init_hmm_bf, control_em = list(restart = list(times = 5)))
#' hmm_4$logLik # -21675.4
#'
#' # Global optimization via StoGO with L-BFGS as final polisher
#' # This can be slow, use parallel computing by adjusting threads argument
#' # (here threads = 1 for portability issues)
#' set.seed(123)
#' hmm_5 <- fit_model(
#'    init_hmm_bf, em_step = FALSE, global_step = TRUE, local_step = TRUE, lb = -50, ub = 50,
#' control_global = list(algorithm = "NLOPT_GD_STOGO", maxeval = 2500, maxtime = 0), threads = 1)
#' hmm_5$logLik # -21675.4
#' 
#' ##############################################################
#' 
#' # Mixture HMM
#'
#' data("biofam3c")
#'
#' ## Building sequence objects
#' marr_seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child_seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left_seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' ## Choosing colors
#' attr(marr_seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(child_seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(left_seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#'
#' ## Starting values for emission probabilities
#' # Cluster 1
#' B1_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.3, 0.6, 0.1, # High probability for married
#'     0.3, 0.3, 0.4), # High probability for divorced
#'   nrow = 4, ncol = 3, byrow = TRUE)
#'
#' B1_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.9, 0.1),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' B1_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9, # High probability for having left home
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' # Cluster 2
#'
#' B2_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.7, 0.2, 0.1),
#'   nrow = 4, ncol = 3, byrow = TRUE)
#'
#' B2_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' B2_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' # Cluster 3
#' B3_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.3, 0.4, 0.3,
#'     0.1, 0.1, 0.8), # High probability for divorced
#'   nrow = 6, ncol = 3, byrow = TRUE)
#'
#' B3_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9),
#'   nrow = 6, ncol = 2, byrow = TRUE)
#'
#'
#' B3_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 6, ncol = 2, byrow = TRUE)
#'
#' # Starting values for transition matrices
#' A1 <- matrix(
#'   c(0.80, 0.16, 0.03, 0.01,
#'     0,    0.90, 0.07, 0.03,
#'     0,    0,    0.90, 0.10,
#'     0,    0,    0,       1),
#'   nrow = 4, ncol = 4, byrow = TRUE)
#'
#' A2 <- matrix(
#'   c(0.80, 0.10, 0.05, 0.03, 0.01, 0.01,
#'     0,    0.70, 0.10, 0.10, 0.05, 0.05,
#'     0,    0,    0.85, 0.01, 0.10, 0.04,
#'     0,    0,    0,    0.90, 0.05, 0.05,
#'     0,    0,    0,    0,    0.90, 0.10,
#'     0,    0,    0,    0,    0,       1),
#'   nrow = 6, ncol = 6, byrow = TRUE)
#'
#' # Starting values for initial state probabilities
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#'
#' # Birth cohort
#' biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
#' biofam3c$covariates$cohort <- factor(
#'   biofam3c$covariates$cohort, labels=c("1909-1935", "1936-1945", "1946-1957"))
#'
#' # Build mixture HMM
#' init_mhmm_bf <- build_mhmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   transition_probs = list(A1, A1, A2),
#'   emission_probs = list(list(B1_marr, B1_child, B1_left),
#'     list(B2_marr, B2_child, B2_left),
#'     list(B3_marr, B3_child, B3_left)),
#'   formula = ~sex + cohort, data = biofam3c$covariates,
#'   channel_names = c("Marriage", "Parenthood", "Residence"))
#'
#'
#' # Fitting the model with different settings
#'
#' # Only EM with default values
#' mhmm_1 <- fit_model(init_mhmm_bf)
#' mhmm_1$logLik # -12713.08
#'
#' # Only L-BFGS
#' mhmm_2 <- fit_model(init_mhmm_bf, em_step = FALSE, local_step = TRUE)
#' mhmm_2$logLik # -12966.51
#'
#' # Use EM with multiple restarts
#' set.seed(123)
#' mhmm_3 <- fit_model(init_mhmm_bf, control_em = list(restart = list(times = 5, transition = FALSE)))
#' mhmm_3$logLik # -12713.08
#' }
#'

fit_model <- function(model, em_step = TRUE, global_step = FALSE, local_step = FALSE,
  control_em = list(), control_global = list(), control_local = list(), lb, ub, threads = 1,
  log_space = FALSE, ...){


  if (!inherits(model, c("hmm", "mhmm"))){
    stop("Argument model must be an object of class 'hmm' or 'mhmm.")
  }

  if (!em_step && !global_step && !local_step) {
    stop("No method chosen for estimation. Choose at least one from em_step, global_step, and local_step.")
  }
  if (threads < 1) {
    stop("Argument threads must be a positive integer.")
  }

  mhmm <- inherits(model, c("mhmm"))

  if (mhmm) {
    df <- attr(model, "df")
    nobs <- attr(model, "nobs")
    original_model <- model
    model <- combine_models(model)
  }

  if (model$n_channels == 1) {
    model$observations <- list(model$observations)
    model$emission_probs <- list(model$emission_probs)
  }

  obsArray <- array(0, c(model$n_sequences, model$length_of_sequences,
    model$n_channels))
  for(i in 1:model$n_channels) {
    obsArray[,,i] <- data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i] > model$n_symbols[i]] <- model$n_symbols[i]
  }
  obsArray <- aperm(obsArray)

  emissionArray<-array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
  for(i in 1:model$n_channels)
    emissionArray[,1:model$n_symbols[i],i]<-model$emission_probs[[i]]

  if (em_step) {
    em.con <- list(print_level = 0, maxeval = 1000, reltol = 1e-12, restart = NULL)
    nmsC <- names(em.con)
    em.con[(namc <- names(control_em))] <- control_em
    if (length(noNms <- namc[!namc %in% nmsC]))
      warning("Unknown names in control_em: ", paste(noNms, collapse = ", "))


    if (!log_space) {
      if (mhmm) {
        resEM <- EMx(model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters,
          em.con$maxeval, em.con$reltol,em.con$print_level, threads)
      } else {
        resEM<-EM(model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, em.con$maxeval, em.con$reltol,em.con$print_level, threads)
      }
    } else {
      if (mhmm) {
        resEM <- log_EMx(model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters,
          em.con$maxeval, em.con$reltol,em.con$print_level, threads)
      } else {
        resEM <- log_EM(model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, em.con$maxeval, em.con$reltol,em.con$print_level, threads)
      }
    }


    if (resEM$error != 0) {
      err_msg <- switch(resEM$error,
        "Scaling factors contain non-finite values.",
        "Backward probabilities contain non-finite values.",
        "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
        "Estimation of coefficients of covariates failed due to singular Hessian.",
        "Estimation of coefficients of covariates failed due to non-finite cluster probabilities.",
        "Non-finite log-likelihood")
      if (!global_step && !local_step &&
          (is.null(control_em$restart$times) || control_em$restart$times == 0)) {
        stop(paste("EM algorithm failed:", err_msg))
      } else warning(paste("EM algorithm failed:", err_msg))
      resEM$logLik <- -Inf
    }


    if (!is.null(em.con$restart)) {
      restart.con <- list(times = 0, print_level = em.con$print_level, maxeval = 100, reltol = 1e-8,
        transition = TRUE, emission = TRUE, sd = 0.25,
        n_optimum = min(control_em$restart$times + 1, 25))
      nmsC <- names(restart.con)
      restart.con[(namc <- names(control_em$restart))] <- control_em$restart
      if (length(noNms <- namc[!namc %in% nmsC]))
        warning("Unknown names in control_em$restart: ", paste(noNms, collapse = ", "))
      restart.con$n_optimum <- min(restart.con$n_optimum, restart.con$times + 1)
    }


    if (!is.null(em.con$restart) && restart.con$times > 0 &&
        (restart.con$transition | restart.con$emission)) {

      opt_restart <- numeric(restart.con$times + 1)
      opt_restart[1] <- resEM$logLik

      if (resEM$error == 0) {
        random_emiss <- resEM$emissionArray
        random_emiss[(random_emiss < 1e-4) & (emissionArray > 0)] <- 1e-4
        for (j in 1:model$n_channels) {
          random_emiss[,1:model$n_symbols[j],j] <-
            random_emiss[,1:model$n_symbols[j],j] / rowSums(random_emiss[,1:model$n_symbols[j],j])
        }
        random_trans <- resEM$transitionMatrix
        random_trans[(random_trans < 1e-4) & (model$transition_probs > 0)] <- 1e-4
        random_trans <- random_trans / rowSums(random_trans)
      } else {
        random_emiss <- emissionArray
        random_trans <- model$transition_probs
      }

      if (restart.con$transition) {
        nz_trans <- (random_trans > 0 & random_trans < 1)
        np_trans <- sum(nz_trans)
        base_trans <- random_trans[nz_trans]
      }
      if (restart.con$emission) {
        nz_emiss <- (random_emiss > 0 & random_emiss < 1)
        np_emiss <- sum(nz_emiss)
        base_emiss <- random_emiss[nz_emiss]
      }

      for (i in 1:restart.con$times) {
        if (restart.con$transition) {
          random_trans[nz_trans] <- abs(base_trans + rnorm(np_trans, sd = restart.con$sd))
          random_trans[(random_trans < 1e-4) & (model$transition_probs > 0)] <- 1e-4
          random_trans <- random_trans / rowSums(random_trans)
        }
        if (restart.con$emission) {
          random_emiss[nz_emiss] <- abs(base_emiss + rnorm(np_emiss, sd = restart.con$sd))
          random_emiss[(random_emiss < 1e-4) & (emissionArray > 0)] <- 1e-4
          for (j in 1:model$n_channels) {
            random_emiss[,1:model$n_symbols[j],j] <-
              random_emiss[,1:model$n_symbols[j],j] / rowSums(random_emiss[,1:model$n_symbols[j],j])
          }
        }
        if (!log_space) {
          if (mhmm) {
            resEMi <- EMx(random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters, restart.con$maxeval,
              restart.con$reltol,restart.con$print_level, threads)
          } else {
            resEMi <- EM(random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, restart.con$maxeval, restart.con$reltol,restart.con$print_level, threads)
          }
        } else {
          if (mhmm) {
            resEMi <- log_EMx(random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters, restart.con$maxeval,
              restart.con$reltol,restart.con$print_level, threads)
          } else {
            resEMi <- log_EM(random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, restart.con$maxeval, restart.con$reltol,restart.con$print_level, threads)
          }
        }


        if (resEMi$error != 0) {
          opt_restart[i + 1] <- -Inf
          err_msg <- switch(resEMi$error,
            "Scaling factors contain non-finite values.",
            "Backward probabilities contain non-finite values.",
            "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
            "Estimation of coefficients of covariates failed due to singular Hessian.",
            "Estimation of coefficients of covariates failed due to non-finite cluster probabilities.",
            "Non-finite log-likelihood")
          warning(paste("EM algorithm failed:", err_msg))
        } else {
          opt_restart[i + 1] <- resEMi$logLik
          if (is.finite(resEMi$logLik) && resEMi$logLik > resEM$logLik) {
            resEM <- resEMi
          }
        }
      }

      opts <- sort(opt_restart, decreasing = TRUE)[1:restart.con$n_optimum]

      if (em.con$reltol < resEM$change) {
        if (!log_space) {
          if (mhmm) {
            resEM <- EMx(resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, resEM$coef, model$X, model$n_states_in_clusters, em.con$maxeval,
              em.con$reltol,em.con$print_level, threads)
          } else {
            resEM <- EM(resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, em.con$maxeval, em.con$reltol,em.con$print_level, threads)
          }

        } else {
          if (mhmm) {
            resEM <- log_EMx(resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, resEM$coef, model$X, model$n_states_in_clusters, em.con$maxeval,
              em.con$reltol,em.con$print_level, threads)
          } else {
            resEM <- log_EM(resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, em.con$maxeval, em.con$reltol, em.con$print_level, threads)
          }
        }
      }

    } else {
      opts <- resEM$logLik
    }

    if (resEM$error == 0) {
      if (resEM$change < -em.con$reltol)
        warning("EM algorithm stopped due to decreasing log-likelihood. ")

      emissionArray <- resEM$emissionArray
      for (i in 1:model$n_channels)
        model$emission_probs[[i]][] <- emissionArray[ , 1:model$n_symbols[i], i]

      if (mhmm && (global_step || local_step)) {
        k <- 0
        for (m in 1:model$n_clusters) {
          original_model$initial_probs[[m]][] <-
            resEM$initialProbs[(k + 1):(k + model$n_states_in_clusters[m])]
          k <- sum(model$n_states_in_clusters[1:m])
        }
      } else {
        model$initial_probs[] <- resEM$initialProbs
      }

      resEM$best_opt_restart <- opts
      model$transition_probs[] <- resEM$transitionMatrix
      model$coefficients[] <- resEM$coefficients
      ll <- resEM$logLik
    } else {
      resEM <- NULL
      ll <- -Inf
    }
  } else resEM <- NULL

  if (global_step || local_step){

    # Largest transition probabilities (for each row)
    x <- which(model$transition_probs > 0, arr.ind = TRUE)
    transNZ <- x[order(x[,1]),]
    maxTM <- cbind(1:model$n_states, max.col(model$transition_probs, ties.method = "first"))
    maxTMvalue <- apply(model$transition_probs,1,max)
    paramTM <- rbind(transNZ,maxTM)
    paramTM <- paramTM[!(duplicated(paramTM) | duplicated(paramTM, fromLast = TRUE)), , drop = FALSE]
    npTM <- nrow(paramTM)
    transNZ <- model$transition_probs > 0
    transNZ[maxTM] <- 0


    emissNZ <- lapply(model$emission_probs, function(i) {
      x <- which(i > 0, arr.ind = TRUE)
      x[order(x[, 1]), ]
    })

    if (model$n_states > 1) {
      maxEM <- lapply(model$emission_probs, function(i){
        cbind(1:model$n_states,max.col(i, ties.method = "first"))
      })
      paramEM <- lapply(1:model$n_channels,function(i) {
        x <- rbind(emissNZ[[i]],maxEM[[i]])
        x[!(duplicated(x)|duplicated(x,fromLast = TRUE)), , drop = FALSE]
      })
      npEM <- sapply(paramEM,nrow)
    } else {
      maxEM <- lapply(model$emission_probs,function(i) max.col(i,ties.method = "first"))
      paramEM <- lapply(1:model$n_channels,function(i) {
        x <- rbind(emissNZ[[i]],c(1,maxEM[[i]]))
        x[!(duplicated(x)|duplicated(x,fromLast = TRUE)),2]
      })
      npEM <- sapply(paramEM, length)
    }
    maxEMvalue <- lapply(1:model$n_channels, function(i)
      apply(model$emission_probs[[i]],1,max))

    emissNZ <- array(0, c(model$n_states,max(model$n_symbols), model$n_channels))
    for (i in 1:model$n_channels) {
      emissNZ[,1:model$n_symbols[i],i]<-model$emission_probs[[i]] > 0
      emissNZ[,1:model$n_symbols[i],i][maxEM[[i]]]<-0
    }

    if (mhmm) {
      maxIP <- maxIPvalue <- npIP <- numeric(original_model$n_clusters)
      paramIP <-  initNZ <- vector("list",original_model$n_clusters)
      for(m in 1:original_model$n_clusters){
        # Index of largest initial probability
        maxIP[m] <- which.max(original_model$initial_probs[[m]])
        # Value of largest initial probability
        maxIPvalue[m] <- original_model$initial_probs[[m]][maxIP[m]]
        # Rest of non-zero probs
        paramIP[[m]] <- setdiff(which(original_model$initial_probs[[m]]>0),maxIP[m])
        npIP[m] <- length(paramIP[[m]])
        initNZ[[m]]<-original_model$initial_probs[[m]] > 0
        initNZ[[m]][maxIP[m]] <- 0
      }
      initNZ<-unlist(initNZ)
      npIPAll <- sum(unlist(npIP))

      npCoef<-length(model$coefficients[,-1])
      model$coefficients[,1] <- 0

      initialvalues<-c(if((npTM+sum(npEM)+npIPAll)>0) log(c(
        if(npTM>0) model$transition_probs[paramTM],
        if(sum(npEM)>0) unlist(sapply(1:model$n_channels,
          function(x) model$emission_probs[[x]][paramEM[[x]]])),
        if(npIPAll>0) unlist(sapply(1:original_model$n_clusters,function(m)
          if(npIP[m]>0) original_model$initial_probs[[m]][paramIP[[m]]]))
      )),
        model$coefficients[,-1]
      )

    } else {
      maxIP <- which.max(model$initial_probs)
      maxIPvalue <- model$initial_probs[maxIP]
      paramIP <- setdiff(which(model$initial_probs > 0), maxIP)

      npIP <- length(paramIP)
      initNZ <- model$initial_probs > 0
      initNZ[maxIP] <- 0

      initialvalues<-c(if((npTM+sum(npEM)+npIP)>0) log(c(
        if(npTM>0) model$transition_probs[paramTM],
        if(sum(npEM)>0) unlist(sapply(1:model$n_channels,
          function(x) model$emission_probs[[x]][paramEM[[x]]])),
        if(npIP>0) model$initial_probs[paramIP]))
      )
    }


    objectivef_mhmm<-function(pars, model, estimate = TRUE){

      if(any(!is.finite(exp(pars))) && estimate)
        return(list(objective = Inf, gradient = rep(-Inf, length(pars))))
      if(npTM>0){
        model$transition_probs[maxTM]<-maxTMvalue
        model$transition_probs[paramTM]<-exp(pars[1:npTM])
        model$transition_probs<-model$transition_probs/rowSums(model$transition_probs)
      }
      if(sum(npEM)>0){
        for(i in 1:model$n_channels){
          emissionArray[,1:model$n_symbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]
          emissionArray[,1:model$n_symbols[i],i][paramEM[[i]]]<-
            exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])
          emissionArray[,1:model$n_symbols[i],i]<-
            emissionArray[,1:model$n_symbols[i],i]/rowSums(emissionArray[,1:model$n_symbols[i],i])
        }
      }
      for(m in 1:original_model$n_clusters){
        if(npIP[m]>0){
          original_model$initial_probs[[m]][maxIP[[m]]] <- maxIPvalue[[m]]
          original_model$initial_probs[[m]][paramIP[[m]]] <- exp(pars[npTM+sum(npEM)+c(0,cumsum(npIP))[m]+
              1:npIP[m]])
          original_model$initial_probs[[m]][] <-
            original_model$initial_probs[[m]]/sum(original_model$initial_probs[[m]])
        }
      }
      model$initial_probs <- unlist(original_model$initial_probs)
      model$coefficients[,-1] <- pars[npTM+sum(npEM)+npIPAll+1:npCoef]

      if (estimate) {
        if (!log_space) {
          if (need_grad) {
            objectivex(model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols,
              model$coefficients, model$X, model$n_states_in_clusters, threads)
          } else {
            -sum(logLikMixHMM(model$transition_probs, emissionArray, model$initial_probs, obsArray,
              model$coefficients, model$X, model$n_states_in_clusters, threads))
          }
        } else {
          if (need_grad) {
            log_objectivex(model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols,
              model$coefficients, model$X, model$n_states_in_clusters, threads)

          } else {
            -sum(log_logLikMixHMM(model$transition_probs, emissionArray, model$initial_probs, obsArray,
              model$coefficients, model$X, model$n_states_in_clusters, threads))
          }
        }

      } else {
        if (sum(npEM) > 0) {
          for (i in 1:model$n_channels) {
            model$emission_probs[[i]][] <- emissionArray[, 1:model$n_symbols[i], i]
          }
        }
        model
      }
    }

    objectivef_hmm <- function(pars, model, estimate = TRUE){

      if(any(!is.finite(exp(pars))) && estimate)
        return(list(objective = Inf, gradient = rep(-Inf, length(pars))))
      if(npTM>0){
        model$transition_probs[maxTM]<-maxTMvalue
        model$transition_probs[paramTM]<-exp(pars[1:npTM])
        model$transition_probs<-model$transition_probs/rowSums(model$transition_probs)
      }
      if(sum(npEM)>0){
        for(i in 1:model$n_channels){
          emissionArray[,1:model$n_symbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]
          emissionArray[,1:model$n_symbols[i],i][paramEM[[i]]]<-
            exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])
          emissionArray[,1:model$n_symbols[i],i]<-
            emissionArray[,1:model$n_symbols[i],i]/rowSums(emissionArray[,1:model$n_symbols[i],i, drop = FALSE])
        }
      }

      if(npIP>0){
        model$initial_probs[maxIP]<-maxIPvalue
        model$initial_probs[paramIP]<-exp(pars[(npTM+sum(npEM)+1):(npTM+sum(npEM)+npIP)])
        model$initial_probs[]<-model$initial_probs/sum(model$initial_probs)
      }


      if (estimate) {
        if (!log_space) {
          if (need_grad) {
            objective(model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols, threads)
          } else {
            -sum(logLikHMM(model$transition_probs, emissionArray,
              model$initial_probs, obsArray, threads))
          }
        } else {
          if (need_grad) {
            log_objective(model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols, threads)
          } else {
            -sum(log_logLikHMM(model$transition_probs, emissionArray,
              model$initial_probs, obsArray, threads))
          }
        }

      } else {
        if (sum(npEM) > 0) {
          for (i in 1:model$n_channels) {
            model$emission_probs[[i]][] <- emissionArray[, 1:model$n_symbols[i], i]
          }
        }
        model
      }
    }

    if (global_step) {

      if (missing(lb)) {
        if (mhmm) {
          lb <- c(rep(-25, length(initialvalues) - npCoef),
            rep(-150 / apply(abs(model$X), 2, max), model$n_clusters - 1))
        } else lb <- -25

      }
      lb <- pmin(lb, 2*initialvalues)
      if (missing(ub)) {
        if (mhmm) {
          ub <- c(rep(25, length(initialvalues) - npCoef),
            rep(150 / apply(abs(model$X), 2, max), model$n_clusters - 1))
        } else ub <- 25
      }
      ub <- pmin(pmax(ub, 2*initialvalues),500)

      if (is.null(control_global$maxeval)) {
        control_global$maxeval <- 10000
      }
      if (is.null(control_global$maxtime)) {
        control_global$maxtime <- 60
      }
      if (is.null(control_global$algorithm)) {
        control_global$algorithm <- "NLOPT_GD_MLSL_LDS"
        if (is.null(control_global$local_opts)) {
          control_global$local_opts <- list(algorithm = "NLOPT_LD_LBFGS",
            ftol_rel = 1e-6, xtol_rel = 1e-4)
        }
      }
      need_grad <- grepl("NLOPT_GD_",control_global$algorithm)

      if (mhmm) {
        globalres <- nloptr(x0 = initialvalues, eval_f = objectivef_mhmm, lb = lb, ub = ub,
          opts = control_global, model = model, estimate = TRUE, ...)
        model <- objectivef_mhmm(globalres$solution, model, FALSE)
      } else {
        globalres <- nloptr(x0 = initialvalues, eval_f = objectivef_hmm, lb = lb, ub = ub,
          opts = control_global, model = model, estimate = TRUE,  ...)
        model <- objectivef_hmm(globalres$solution, model, FALSE)
      }

      if (globalres$status < 0) {
        warning(paste("Global optimization terminated:", globalres$message))
      }
      initialvalues <- globalres$solution
      ll <- -globalres$objective

    } else globalres <- NULL

    if (local_step) {
      if ( is.null(control_local$maxeval)) {
        control_local$maxeval <- 10000
      }
      if (is.null(control_local$algorithm)) {
        control_local$algorithm <- "NLOPT_LD_LBFGS"
      }
      if (is.null(control_local$xtol_rel)) {
        control_local$xtol_rel <- 1e-8
      }
      if (is.null(control_local$ftol_rel)) {
        control_local$ftol_rel <- 1e-10
      }

      need_grad <- grepl("NLOPT_LD_",control_local$algorithm)
      if (mhmm) {
        localres <- nloptr(x0 = initialvalues, eval_f = objectivef_mhmm,
          opts = control_local, model = model, estimate = TRUE, ...)
        model <- objectivef_mhmm(localres$solution, model, FALSE)

      } else {
        localres <- nloptr(x0 = initialvalues, eval_f = objectivef_hmm,
          opts = control_local, model = model, estimate = TRUE, ...)
        model <- objectivef_hmm(localres$solution, model, FALSE)
      }
      if (localres$status < 0) {
        warning(paste("Local optimization terminated:", localres$message))
      }
      ll <- -localres$objective
    } else localres <- NULL

    if (mhmm) {
      rownames(model$coefficients) <- colnames(model$X)
      colnames(model$coefficients) <- model$cluster_names
    }
  } else globalres <- localres <- NULL

  if (model$n_channels == 1) {
    model$observations <- model$observations[[1]]
    model$emission_probs <- model$emission_probs[[1]]
  }

  if (mhmm) {
    model <- spread_models(model)
    attr(model, "df") <- df
    attr(model, "nobs") <- nobs
    attr(model, "type") <- attr(original_model, "type")

    for (i in 1:model$n_clusters) {
      dimnames(model$transition_probs[[i]]) <- dimnames(original_model$transition_probs[[i]])
      for (j in 1:model$n_channels) {
        dimnames(model$emission_probs[[i]][[j]]) <- dimnames(original_model$emission_probs[[i]][[j]])
      }
    }
  }

  suppressWarnings(try(model <- trim_model(model, verbose = FALSE), silent = TRUE))
  list(model = model,
    logLik = ll, em_results = resEM[c("logLik", "iterations", "change", "best_opt_restart")],
    global_results = globalres, local_results = localres, call = match.call())

}
