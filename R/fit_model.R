#' Estimate Parameters of (Mixture) Hidden Markov Models and Their Restricted
#' Variants
#'
#' Function `fit_model` estimates the parameters of mixture hidden
#' Markov models and its restricted variants using maximum likelihood.
#' Initial values for estimation are taken from the corresponding components
#' of the model with preservation of original zero probabilities.
#'
#' @export
#' @param model An object of class `hmm` or `mhmm`.
#' @param em_step Logical. Whether or not to use the EM algorithm at the start
#'   of the parameter estimation. The default is `TRUE`.
#' @param global_step Logical. Whether or not to use global optimization via
#'   [nloptr::nloptr()] (possibly after the EM step). The default is `FALSE`.
#' @param local_step Logical. Whether or not to use local optimization via
#'   [nloptr::nloptr()] (possibly after the EM and/or global steps). The 
#'   default is `FALSE`.
#' @param control_em Optional list of control parameters for the EM algorithm.
#' Possible arguments are 
#' * `maxeval`\cr The maximum number of iterations, the default is 1000. Note that iteration 
#'   counter starts with -1 so with `maxeval = 1` you get already two iterations.
#'   This is for backward compatibility reasons.
#' * `print_level`\cr The level of printing. Possible values are 0
#'   (prints nothing), 1 (prints information at the start and the end of the 
#'   algorithm), 2 (prints at every iteration), and for mixture models 3 
#'   (print also during optimization of coefficients).
#' * `reltol`\cr Relative tolerance for convergence defined as 
#'   \eqn{(logLik_new - logLik_old)/(abs(logLik_old) + 0.1)}. The default is 1e-10.
#' * `restart`\cr A list containing options for possible EM restarts with the
#'   following components:
#'   * `times`\cr Number of restarts of the EM algorithm using random initial 
#'     values. The default is `0`, i.e. no restarts.
#'   * `transition`\cr Logical. Should the original transition probabilities be 
#'     varied? The default is `TRUE`.
#'   * `emission`\cr Logical. Should the original emission probabilities be 
#'     varied? The default is `TRUE`.
#'   * `sd`\cr Standard deviation for [stats::rnorm()] used in randomization. The 
#'     default is 0.25.
#'   * `maxeval`\cr Maximum number of iterations, the default is 
#'     `control_em$maxeval`
#'   * print_level\cr Level of printing in restarted EM steps. The default is 
#'     `control_em$print_level`.
#'   * `reltol`\cr Relative tolerance for convergence at restarted EM steps. 
#'     The default is `control_em$reltol`. If the relative change of the final 
#'     model of the restart phase is larger than the tolerance for the original 
#'     EM phase, the final model is re-estimated with the original `reltol`
#'     and `maxeval` at the end of the EM step.
#'   * `n_optimum`\cr Save the log-likelihood values of the `n_optimum` best
#'     models (from all estimated models including the the first EM run.).
#'     The default is `min(times + 1, 25)`.
#'   * `use_original`\cr If `TRUE`, use the initial values of the input model as 
#'     starting points for the permutations. Otherwise permute the results of 
#'     the first EM run.
#' @param control_global Optional list of additional arguments for 
#' [nloptr::nloptr()] argument `opts`. The default values are
#'  * `algorithm`\cr `"NLOPT_GD_MLSL_LDS"`
#'  * `local_opts`\cr `list(algorithm = "NLOPT_LD_LBFGS", ftol_rel = 1e-6, xtol_rel = 1e-4)`
#'  * `maxeval`\cr `10000` (maximum number of iterations in global 
#'    optimization algorithm.)
#'  * `maxtime`\cr `60` (maximum time for global optimization. Set to 0 for 
#'    unlimited time.)
#' @param lb,ub Lower and upper bounds for parameters in Softmax 
#' parameterization.The default interval is 
#' `c(pmin(-25, 2*initialvalues), pmax(25, 2*initialvalues))`, except for gamma 
#' coefficients,where the scale of covariates is taken into account.
#' Note that it might still be a good idea to scale covariates around unit 
#' scale. Bounds are used only in the global optimization step.
#' @param control_local Optional list of additional arguments for
#'   [nloptr::nloptr()] argument `opts`. The default values are
#'  * `algorithm`\cr `"NLOPT_LD_LBFGS"`
#'  * `ftol_rel`\cr `1e-10`
#'  * `xtol_rel`\cr `1e-8`
#'  * `maxeval`\cr `10000` (maximum number of iterations)
#' @param threads Number of threads to use in parallel computing. The default 
#' is `1`.
#' @param log_space Make computations using log-space instead of scaling for 
#' greater numerical stability at a cost of decreased computational performance. 
#' The default is `FALSE` (but was `TRUE` in 2.0.0 version of the package).
#' @param constraints Integer vector defining equality constraints for emission 
#' distributions. Not supported for EM algorithm. See details.
#' @param fixed_inits Can be used to fix some of the probabilities to their 
#' initial values.Should have same structure as `model$initial_probs`, where 
#' each element is either `TRUE` (fixed) or `FALSE` (to be estimated). Note that 
#' zero probabilities are always fixed to 0. Not supported for EM algorithm. 
#' See details.
#' @param fixed_emissions  Can be used to fix some of the probabilities to 
#' their initial values. Should have same structure as `model$emission_probs`, 
#' where each element is either `TRUE` (fixed) or `FALSE` (to be estimated). Note 
#' that zero probabilities are always fixed to 0. Not supported for EM 
#' algorithm. See details.
#' @param fixed_transitions  Can be used to fix some of the probabilities to 
#' their initial values. Should have same structure as `model$transition_probs`, 
#' where each element is either `TRUE` (fixed) or `FALSE` (to be estimated). 
#' Note that zero probabilities are always fixed to 0. Not supported for EM 
#' algorithm. See details.
#' @param ... Additional arguments to [nloptr::nloptr()].
#' @return
#' * logLik\cr Log-likelihood of the estimated model.
#' * em_results\cr Results after the EM step: log-likelihood (`logLik`), 
#'   number of iterations (`iterations`), relative change in log-likelihoods 
#'   between the last two iterations (`change`), and the log-likelihoods of 
#'   the `n_optimum` best models after the EM step (`best_opt_restart`).
#' * global_results\cr Results after the global step.
#' * local_results\cr Results after the local step.
#' * call\cr The matched function call.
#' @seealso [build_hmm()],  [build_mhmm()],
#' [build_mm()],  [build_mmm()], and  [build_lcm()]
#'   for constructing different types of models; [summary.mhmm()]
#'   for a summary of a MHMM; [separate_mhmm()] for reorganizing a MHMM into
#'   a list of separate hidden Markov models; and [plot.hmm()] and [plot.mhmm()]
#'   for plotting model objects.
#' @details The fitting function provides three estimation steps: 1) EM algorithm,
#'   2) global optimization, and 3) local optimization. The user can call for one method
#'   or any combination of these steps, but should note that they are preformed in the
#'   above-mentioned order. The results from a former step are used as starting values
#'   in a latter, except for some of global optimization algorithms (such as MLSL and StoGO)
#'   which only use initial values for setting up the boundaries for the optimization.
#'
#'   It is possible to rerun the EM algorithm automatically using random starting
#'   values based on the first run of EM. Number of restarts is defined by
#'   the `restart` argument in `control_em`. As the EM algorithm is
#'   relatively fast, this method might be preferred option compared to the proper
#'   global optimization strategy of step 2.
#'
#'   The default global optimization method (triggered via `global_step = TRUE`) is
#'   the multilevel single-linkage method (MLSL) with the LDS modification (`NLOPT_GD_MLSL_LDS` as
#'   `algorithm` in `control_global`), with L-BFGS as the local optimizer.
#'   The MLSL method draws random starting points and performs a local optimization
#'   from each. The LDS modification uses low-discrepancy sequences instead of
#'   pseudo-random numbers as starting points and should improve the convergence rate.
#'   In order to reduce the computation time spent on non-global optima, the
#'   convergence tolerance of the local optimizer is set relatively large. At step 3,
#'   a local optimization (L-BFGS by default) is run with a lower tolerance to find the
#'   optimum with high precision.
#'
#'   There are some theoretical guarantees that the MLSL method used as the default
#'   optimizer in step 2 should find all local optima in a finite number of local
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
#'   Any algorithm available in the `nloptr` function can be used for the global and
#'   local steps.
#'
#'   Equality constraints for emission distributions can be defined using the argument
#'   `constraints`. This should be a vector with length equal to the number of states,
#'   with numbers starting from 1 and increasing for each unique row of the emission probability matrix.
#'   For example in case of five states with emissions of first and third states being equal,
#'   `constraints = c(1, 2, 1, 3, 4)`. Similarly, some of the model parameters can be fixed to their
#'   initial values by using arguments `fixed_inits`, `fixed_emissions`,
#'   and `fixed_transitions`, where the structure of the arguments should be
#'   same as the corresponding model components, so that TRUE value means that
#'   the parameter should be fixed and FALSE otherwise (it is still treated as fixed if it
#'   is zero though). For both types of constrains, only numerical optimisation
#'   (local or global) is available, and currently the gradients are computed numerically
#'   (if needed) in these cases.
#'
#'   In a case where the is no transitions from one state to anywhere (even to
#'   itself), the state is defined as absorbing in a way that probability of
#'   staying in this state is fixed to 1. See also `build_mm` function.
#'
#' @references Helske S. and Helske J. (2019). Mixture Hidden Markov Models for Sequence Data: The seqHMM Package in R,
#' Journal of Statistical Software, 88(3), 1-32. doi:10.18637/jss.v088.i03
#'
#' @examples
#'
#' # Hidden Markov model for mvad data
#'
#' data("mvad", package = "TraMineR")
#'
#' mvad_alphabet <-
#'   c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad_labels <- c(
#'   "employment", "further education", "higher education",
#'   "joblessness", "school", "training"
#' )
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 15:86,
#'   alphabet = mvad_alphabet,
#'   states = mvad_scodes, labels = mvad_labels, xtstep = 6,
#'   cpal = colorpalette[[6]]
#' )
#'
#' # Starting values for the emission matrix
#' emiss <- matrix(
#'   c(
#'     0.05, 0.05, 0.05, 0.05, 0.75, 0.05, # SC
#'     0.05, 0.75, 0.05, 0.05, 0.05, 0.05, # FE
#'     0.05, 0.05, 0.05, 0.4, 0.05, 0.4, # JL, TR
#'     0.05, 0.05, 0.75, 0.05, 0.05, 0.05, # HE
#'     0.75, 0.05, 0.05, 0.05, 0.05, 0.05
#'   ), # EM
#'   nrow = 5, ncol = 6, byrow = TRUE
#' )
#'
#' # Starting values for the transition matrix
#' trans <- matrix(0.025, 5, 5)
#' diag(trans) <- 0.9
#'
#' # Starting values for initial state probabilities
#' initial_probs <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#'
#' # Building a hidden Markov model
#' init_hmm_mvad <- build_hmm(
#'   observations = mvad_seq,
#'   transition_probs = trans, emission_probs = emiss,
#'   initial_probs = initial_probs
#' )
#'
#' \dontrun{
#' set.seed(21)
#' fit_hmm_mvad <- fit_model(init_hmm_mvad, control_em = list(restart = list(times = 50)))
#' hmm_mvad <- fit_hmm_mvad$model
#' }
#'
#' # save time, load the previously estimated model
#' data("hmm_mvad")
#'
#' # Markov model
#' # Note: build_mm estimates model parameters from observations,
#' # no need for estimating with fit_model unless there are missing observations
#'
#' mm_mvad <- build_mm(observations = mvad_seq)
#'
#' # Comparing likelihoods, MM fits better
#' logLik(hmm_mvad)
#' logLik(mm_mvad)
#'
#' \dontrun{
#' require("igraph") # for layout_in_circle
#'
#' plot(mm_mvad,
#'   layout = layout_in_circle, legend.prop = 0.3,
#'   edge.curved = 0.3, edge.label = NA,
#'   vertex.label.pos = c(0, 0, pi, pi, pi, 0)
#' )
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
#' marr_seq <- seqdef(biofam3c$married,
#'   start = 15,
#'   alphabet = c("single", "married", "divorced"),
#'   cpal = c("violetred2", "darkgoldenrod2", "darkmagenta")
#' )
#' child_seq <- seqdef(biofam3c$children,
#'   start = 15,
#'   alphabet = c("childless", "children"),
#'   cpal = c("darkseagreen1", "coral3")
#' )
#' left_seq <- seqdef(biofam3c$left,
#'   start = 15,
#'   alphabet = c("with parents", "left home"),
#'   cpal = c("lightblue", "red3")
#' )
#'
#' # Starting values for emission matrices
#'
#' emiss_marr <- matrix(NA, nrow = 3, ncol = 3)
#' emiss_marr[1, ] <- seqstatf(marr_seq[, 1:5])[, 2] + 1
#' emiss_marr[2, ] <- seqstatf(marr_seq[, 6:10])[, 2] + 1
#' emiss_marr[3, ] <- seqstatf(marr_seq[, 11:16])[, 2] + 1
#' emiss_marr <- emiss_marr / rowSums(emiss_marr)
#'
#' emiss_child <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_child[1, ] <- seqstatf(child_seq[, 1:5])[, 2] + 1
#' emiss_child[2, ] <- seqstatf(child_seq[, 6:10])[, 2] + 1
#' emiss_child[3, ] <- seqstatf(child_seq[, 11:16])[, 2] + 1
#' emiss_child <- emiss_child / rowSums(emiss_child)
#'
#' emiss_left <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_left[1, ] <- seqstatf(left_seq[, 1:5])[, 2] + 1
#' emiss_left[2, ] <- seqstatf(left_seq[, 6:10])[, 2] + 1
#' emiss_left[3, ] <- seqstatf(left_seq[, 11:16])[, 2] + 1
#' emiss_left <- emiss_left / rowSums(emiss_left)
#'
#' # Starting values for transition matrix
#' trans <- matrix(c(
#'   0.9, 0.07, 0.03,
#'   0, 0.9, 0.1,
#'   0, 0, 1
#' ), nrow = 3, ncol = 3, byrow = TRUE)
#'
#' # Starting values for initial state probabilities
#' inits <- c(0.9, 0.09, 0.01)
#'
#' # Building hidden Markov model with initial parameter values
#' init_hmm_bf <- build_hmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   transition_probs = trans,
#'   emission_probs = list(emiss_marr, emiss_child, emiss_left),
#'   initial_probs = inits
#' )
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
#'   init_hmm_bf,
#'   em_step = FALSE, global_step = TRUE, local_step = TRUE,
#'   control_global = list(maxeval = 5000, maxtime = 0), threads = 1
#' )
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
#'   init_hmm_bf,
#'   em_step = FALSE, global_step = TRUE, local_step = TRUE,
#'   lb = -50, ub = 50, control_global = list(
#'     algorithm = "NLOPT_GD_STOGO",
#'     maxeval = 2500, maxtime = 0
#'   ), threads = 1
#' )
#' hmm_5$logLik # -21675.4
#'
#' ##############################################################
#'
#' # Mixture HMM
#'
#' data("biofam3c")
#'
#' # Building sequence objects
#' marr_seq <- seqdef(biofam3c$married,
#'   start = 15,
#'   alphabet = c("single", "married", "divorced"),
#'   cpal = c("violetred2", "darkgoldenrod2", "darkmagenta")
#' )
#' child_seq <- seqdef(biofam3c$children,
#'   start = 15,
#'   alphabet = c("childless", "children"),
#'   cpal = c("darkseagreen1", "coral3")
#' )
#' left_seq <- seqdef(biofam3c$left,
#'   start = 15,
#'   alphabet = c("with parents", "left home"),
#'   cpal = c("lightblue", "red3")
#' )
#'
#' ## Starting values for emission probabilities
#' # Cluster 1
#' B1_marr <- matrix(
#'   c(
#'     0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.3, 0.6, 0.1, # High probability for married
#'     0.3, 0.3, 0.4
#'   ), # High probability for divorced
#'   nrow = 4, ncol = 3, byrow = TRUE
#' )
#'
#' B1_child <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.9, 0.1
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' B1_left <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9, # High probability for having left home
#'     0.1, 0.9,
#'     0.1, 0.9
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' # Cluster 2
#'
#' B2_marr <- matrix(
#'   c(
#'     0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.7, 0.2, 0.1
#'   ),
#'   nrow = 4, ncol = 3, byrow = TRUE
#' )
#'
#' B2_child <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.1, 0.9
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' B2_left <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.1, 0.9
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' # Cluster 3
#' B3_marr <- matrix(
#'   c(
#'     0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.3, 0.4, 0.3,
#'     0.1, 0.1, 0.8
#'   ), # High probability for divorced
#'   nrow = 6, ncol = 3, byrow = TRUE
#' )
#'
#' B3_child <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9
#'   ),
#'   nrow = 6, ncol = 2, byrow = TRUE
#' )
#'
#'
#' B3_left <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9,
#'     0.1, 0.9
#'   ),
#'   nrow = 6, ncol = 2, byrow = TRUE
#' )
#'
#' # Starting values for transition matrices
#' A1 <- matrix(
#'   c(
#'     0.80, 0.16, 0.03, 0.01,
#'     0, 0.90, 0.07, 0.03,
#'     0, 0, 0.90, 0.10,
#'     0, 0, 0, 1
#'   ),
#'   nrow = 4, ncol = 4, byrow = TRUE
#' )
#'
#' A2 <- matrix(
#'   c(
#'     0.80, 0.10, 0.05, 0.03, 0.01, 0.01,
#'     0, 0.70, 0.10, 0.10, 0.05, 0.05,
#'     0, 0, 0.85, 0.01, 0.10, 0.04,
#'     0, 0, 0, 0.90, 0.05, 0.05,
#'     0, 0, 0, 0, 0.90, 0.10,
#'     0, 0, 0, 0, 0, 1
#'   ),
#'   nrow = 6, ncol = 6, byrow = TRUE
#' )
#'
#' # Starting values for initial state probabilities
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#'
#' # Birth cohort
#' biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
#' biofam3c$covariates$cohort <- factor(
#'   biofam3c$covariates$cohort,
#'   labels = c("1909-1935", "1936-1945", "1946-1957")
#' )
#'
#' # Build mixture HMM
#' init_mhmm_bf <- build_mhmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   transition_probs = list(A1, A1, A2),
#'   emission_probs = list(
#'     list(B1_marr, B1_child, B1_left),
#'     list(B2_marr, B2_child, B2_left),
#'     list(B3_marr, B3_child, B3_left)
#'   ),
#'   formula = ~ sex + cohort, data = biofam3c$covariates,
#'   channel_names = c("Marriage", "Parenthood", "Residence")
#' )
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
#' # Left-to-right HMM with equality constraint:
#'
#' set.seed(1)
#'
#' # Transition matrix
#' # Either stay or move to next state
#' A <- diag(c(0.9, 0.95, 0.95, 1))
#' A[1, 2] <- 0.1
#' A[2, 3] <- 0.05
#' A[3, 4] <- 0.05
#'
#' # Emission matrix, rows 1 and 3 equal
#' B <- rbind(
#'   c(0.4, 0.2, 0.3, 0.1),
#'   c(0.1, 0.5, 0.1, 0.3),
#'   c(0.4, 0.2, 0.3, 0.1),
#'   c(0, 0.2, 0.4, 0.4)
#' )
#'
#' # Start from first state
#' init <- c(1, 0, 0, 0)
#'
#' # Simulate sequences
#' sim <- simulate_hmm(
#'   n_sequences = 100,
#'   sequence_length = 20, init, A, B
#' )
#'
#' # initial model, use true values as inits for faster estimation here
#' model <- build_hmm(sim$observations, init = init, trans = A, emiss = B)
#'
#' # estimate the model subject to constraints:
#' # First and third row of emission matrix are equal (see details)
#' fit <- fit_model(model,
#'   constraints = c(1, 2, 1, 3),
#'   em_step = FALSE, local_step = TRUE
#' )
#' fit$model
#'
#' ## Fix some emissions:
#'
#' fixB <- matrix(FALSE, 4, 4)
#' fixB[2, 1] <- fixB[1, 3] <- TRUE # these are fixed to their initial values
#' fit <- fit_model(model,
#'   fixed_emissions = fixB,
#'   em_step = FALSE, local_step = TRUE
#' )
#' fit$model$emission_probs
fit_model <- function(
    model, em_step = TRUE, global_step = FALSE, local_step = FALSE,
    control_em = list(), control_global = list(), control_local = list(), lb, ub, threads = 1,
    log_space = FALSE, constraints = NULL, fixed_inits = NULL,
    fixed_emissions = NULL, fixed_transitions = NULL, ...) {
  
  stopifnot_(
    checkmate::test_int(x = threads, lower = 1L), 
    "Argument {.arg threads} must be a single positive integer."
  )
  # check that the initial model is ok
  stopifnot_(
    is.finite(logLik(model)),
    "Initial log-likelihood of the input model is nonfinite. Please adjust the 
    initial values."
  )
  stopifnot_(
    inherits(model, c("hmm", "mhmm")),
    "{.arg model} {.cls hmm} or {.cls mhmm} object."
  )
  stopifnot_(
    em_step || global_step || local_step,
    "No method chosen for estimation. Choose at least one from `em_step`, 
          `global_step`, and `local_step`."
  )
  stopifnot_(
    is.null(constraints) || !em_step,
    "EM algorithm does not support constraints. Use local and/or global 
    steps only."
  )
  
  fixed <- !is.null(fixed_inits) | !is.null(fixed_emissions) | !is.null(fixed_transitions)
  stopifnot_(
    !fixed || !em_step,
    "EM algorithm does not support fixed probabilities (except zeros). Use 
    local and/or global steps only."
  )
  mhmm <- inherits(model, c("mhmm"))
  if (mhmm) {
    original_model <- model
    model <- .combine_models(model)
  }
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)
  if (model$n_channels == 1) {
    model$emission_probs <- list(model$emission_probs)
  }
  if (em_step) {
    em.con <- list(print_level = 0, maxeval = 1000, reltol = 1e-10, restart = NULL)
    nmsC <- names(em.con)
    em.con[(namc <- names(control_em))] <- control_em
    if (length(noNms <- namc[!namc %in% nmsC])) {
      warning_("Unknown names in {.arg control_em}: ", 
               paste(noNms, collapse = ", "))
    }
    if (!log_space) {
      if (mhmm) {
        resEM <- EMx(
          model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters,
          em.con$maxeval, em.con$reltol, em.con$print_level, threads
        )
      } else {
        resEM <- EM(
          model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, em.con$maxeval, em.con$reltol, em.con$print_level, threads
        )
      }
    } else {
      if (mhmm) {
        resEM <- log_EMx(
          model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters,
          em.con$maxeval, em.con$reltol, em.con$print_level, threads
        )
      } else {
        resEM <- log_EM(
          model$transition_probs, emissionArray, model$initial_probs, obsArray,
          model$n_symbols, em.con$maxeval, em.con$reltol, em.con$print_level, threads
        )
      }
    }
    
    
    if (resEM$error != 0) {
      err_msg <- switch(resEM$error,
                        "Scaling factors contain non-finite values.",
                        "Backward probabilities contain non-finite values.",
                        "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
                        "Estimation of gamma coefficients failed due to singular Hessian.",
                        "Estimation of gamma coefficients failed due to non-finite cluster probabilities.",
                        "Non-finite log-likelihood"
      )
      if (!global_step && !local_step &&
          (is.null(control_em$restart$times) || control_em$restart$times == 0)) {
        stop_(paste("EM algorithm failed:", err_msg))
      } else {
        warning_(paste("EM algorithm failed:", err_msg))
      }
      resEM$logLik <- -Inf
    }
    
    
    if (!is.null(em.con$restart)) {
      restart.con <- list(
        times = 0, print_level = em.con$print_level,
        maxeval = em.con$maxeval, reltol = em.con$reltol,
        transition = TRUE, emission = TRUE, coef = FALSE, sd = 0.25,
        n_optimum = min(control_em$restart$times + 1, 25), use_original = TRUE
      )
      nmsC <- names(restart.con)
      restart.con[(namc <- names(control_em$restart))] <- control_em$restart
      if (length(noNms <- namc[!namc %in% nmsC])) {
        warning_("Unknown names in {.arg control_em$restart}: 
                 {paste(noNms, collapse = ", ")}")
      }
      restart.con$n_optimum <- min(restart.con$n_optimum, restart.con$times + 1)
    }
    if (!is.null(em.con$restart) && restart.con$times > 0 &&
        (restart.con$transition | restart.con$emission)) {
      opt_restart <- numeric(restart.con$times + 1)
      opt_restart[1] <- resEM$logLik
      
      if (restart.con$use_original) {
        random_emiss <- emissionArray
        random_trans <- model$transition_probs
        random_coef <- model$coefficients
      } else {
        random_emiss <- resEM$emission
        random_trans <- resEM$transition
        random_coef <- resEM$coef
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
      
      if (restart.con$coef) {
        base_coef <- random_coef
        coef_scale <- max(abs(base_coef))
        n_coef <- length(base_coef)
      }
      
      for (i in 1:restart.con$times) {
        if (restart.con$transition) {
          random_trans[nz_trans] <- abs(base_trans + stats::rnorm(np_trans, sd = restart.con$sd))
          random_trans[(random_trans < 1e-4) & (model$transition_probs > 0)] <- 1e-4
          random_trans <- random_trans / rowSums(random_trans)
        }
        if (restart.con$emission) {
          random_emiss[nz_emiss] <- abs(base_emiss + stats::rnorm(np_emiss, sd = restart.con$sd))
          random_emiss[(random_emiss < 1e-4) & (emissionArray > 0)] <- 1e-4
          for (j in 1:model$n_channels) {
            random_emiss[, 1:model$n_symbols[j], j] <-
              random_emiss[, 1:model$n_symbols[j], j] / rowSums(random_emiss[, 1:model$n_symbols[j], j])
          }
        }
        if (restart.con$coef) {
          random_coef[] <- stats::rnorm(n_coef, base_coef, sd = coef_scale)
        }
        
        if (!log_space) {
          if (mhmm) {
            resEMi <- EMx(
              random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, random_coef, model$X, model$n_states_in_clusters, restart.con$maxeval,
              restart.con$reltol, restart.con$print_level, threads
            )
          } else {
            resEMi <- EM(
              random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, restart.con$maxeval, restart.con$reltol, restart.con$print_level, threads
            )
          }
        } else {
          if (mhmm) {
            resEMi <- log_EMx(
              random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, random_coef, model$X, model$n_states_in_clusters, restart.con$maxeval,
              restart.con$reltol, restart.con$print_level, threads
            )
          } else {
            resEMi <- log_EM(
              random_trans, random_emiss, model$initial_probs, obsArray,
              model$n_symbols, restart.con$maxeval, restart.con$reltol, restart.con$print_level, threads
            )
          }
        }
        
        
        if (resEMi$error != 0) {
          opt_restart[i + 1] <- -Inf
          err_msg <- switch(resEMi$error,
                            "Scaling factors contain non-finite values.",
                            "Backward probabilities contain non-finite values.",
                            "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
                            "Estimation of gamma coefficients failed due to singular Hessian.",
                            "Estimation of gamma coefficients failed due to non-finite cluster probabilities.",
                            "Non-finite log-likelihood"
          )
          warning_(paste("EM algorithm failed:", err_msg))
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
            resEM <- EMx(
              resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, resEM$coef, model$X, model$n_states_in_clusters, em.con$maxeval,
              em.con$reltol, em.con$print_level, threads
            )
          } else {
            resEM <- EM(
              resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, em.con$maxeval, em.con$reltol, em.con$print_level, threads
            )
          }
        } else {
          if (mhmm) {
            resEM <- log_EMx(
              resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, resEM$coef, model$X, model$n_states_in_clusters, em.con$maxeval,
              em.con$reltol, em.con$print_level, threads
            )
          } else {
            resEM <- log_EM(
              resEM$transitionMatrix, resEM$emissionArray, resEM$initialProbs, obsArray,
              model$n_symbols, em.con$maxeval, em.con$reltol, em.con$print_level, threads
            )
          }
        }
      }
    } else {
      opts <- resEM$logLik
    }
    
    if (resEM$error == 0) {
      if (resEM$change < -em.con$reltol) {
        warning_("EM algorithm stopped due to decreasing log-likelihood.")
      }
      
      emissionArray <- resEM$emissionArray
      for (i in 1:model$n_channels) {
        model$emission_probs[[i]][] <- emissionArray[, 1:model$n_symbols[i], i]
      }
      
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
  } else {
    resEM <- NULL
  }
  
  if (global_step || local_step) {
    if (is.null(fixed_transitions)) {
      transNZ <- model$transition_probs > 0
    } else {
      if (mhmm) {
        fixed_transitions <- as.matrix(Matrix::.bdiag(fixed_transitions))
      }
      mode(fixed_transitions) <- "logical"
      transNZ <- model$transition_probs > 0 & !fixed_transitions
    }
    # single value per row => nothing to estimate
    transNZ[rowSums(transNZ) == 1, ] <- 0
    # Largest transition probabilities (for each row)
    maxTM <- cbind(1:model$n_states, apply(model$transition_probs, 1, which.max))
    maxTMvalue <- model$transition_probs[maxTM]
    transNZ[maxTM] <- 0
    paramTM <- which(transNZ == 1, arr.ind = TRUE)
    paramTM <- paramTM[order(paramTM[, 1]), ]
    npTM <- nrow(paramTM)
    n_states <- model$n_states
    original_emiss <- model$emission_probs
    original_trans <- model$transition_probs
    if (!is.null(constraints)) {
      stopifnot_(
        identical(length(constraints), model$n_states),
        "Length of {.arg constraints} is not equal to the number of states."
      )
      uniqs <- !duplicated(constraints)
      model$emission_probs <- lapply(model$emission_probs, \(x) {
        x[uniqs, ]
      })
      n_states <- sum(uniqs)
    }
    if (is.null(fixed_emissions)) {
      emissNZ <- lapply(model$emission_probs, \(i) {
        i > 0
      })
    } else {
      if (mhmm) {
        if (model$n_channels > 1) {
          fixed_emissions <- lapply(1:model$n_channels, \(i) {
            do.call(rbind, sapply(fixed_emissions, "[", i))
          })
        } else {
          fixed_emissions <- list(do.call(rbind, fixed_emissions))
        }
      } else {
        if (model$n_channels == 1) {
          fixed_emissions <- list(fixed_emissions)
        }
      }
      for (i in 1:model$n_channels) mode(fixed_emissions[[i]]) <- "logical"
      if (!is.null(constraints)) {
        for (i in 1:model$n_channels) fixed_emissions[[i]] <- fixed_emissions[[i]][uniqs, ]
      }
      emissNZ <- lapply(1:model$n_channels, \(i) {
        model$emission_probs[[i]] > 0 & !fixed_emissions[[i]]
      })
    }
    maxEM <- lapply(1:model$n_channels, \(i) {
      cbind(1:n_states, apply(model$emission_probs[[i]], 1, which.max))
    })
    maxEMvalue <- lapply(1:model$n_channels, \(i) {
      apply(model$emission_probs[[i]], 1, max)
    })
    emissNZ <- array(0, c(n_states, max(model$n_symbols), model$n_channels))
    npEM <- numeric(model$n_channels)
    paramEM <- vector("list", model$n_channels)
    if (is.null(fixed_emissions)) {
      for (i in 1:model$n_channels) {
        emissNZ[, 1:model$n_symbols[i], i] <- model$emission_probs[[i]] > 0
        emissNZ[, 1:model$n_symbols[i], i][maxEM[[i]]] <- 0
        paramEM[[i]] <- which(emissNZ[, 1:model$n_symbols[i], i] == 1, arr.ind = TRUE)
        paramEM[[i]] <- paramEM[[i]][order(paramEM[[i]][, 1]), , drop = FALSE]
        npEM[i] <- nrow(paramEM[[i]])
      }
    } else {
      for (i in 1:model$n_channels) {
        emissNZ[, 1:model$n_symbols[i], i] <- model$emission_probs[[i]] > 0 & !fixed_emissions[[i]]
        emissNZ[, 1:model$n_symbols[i], i][maxEM[[i]]] <- 0
        paramEM[[i]] <- which(emissNZ[, 1:model$n_symbols[i], i] == 1, arr.ind = TRUE)
        paramEM[[i]] <- paramEM[[i]][order(paramEM[[i]][, 1]), , drop = FALSE]
        npEM[i] <- nrow(paramEM[[i]])
      }
    }
    
    if (mhmm) {
      maxIP <- maxIPvalue <- npIP <- numeric(original_model$n_clusters)
      paramIP <- initNZ <- vector("list", original_model$n_clusters)
      for (m in 1:original_model$n_clusters) {
        # Index of largest initial probability
        maxIP[m] <- which.max(original_model$initial_probs[[m]])
        # Value of largest initial probability
        maxIPvalue[m] <- original_model$initial_probs[[m]][maxIP[m]]
        # Rest of non-zero probs
        if (is.null(fixed_inits)) {
          initNZ[[m]] <- original_model$initial_probs[[m]] > 0
        } else {
          initNZ[[m]] <- original_model$initial_probs[[m]] > 0 & !fixed_inits[[m]]
        }
        initNZ[[m]][maxIP[m]] <- 0
        paramIP[[m]] <- which(initNZ[[m]] == 1)
        npIP[m] <- length(paramIP[[m]])
      }
      initNZ <- unlist(initNZ)
      npIPAll <- sum(unlist(npIP))
      
      npCoef <- length(model$coefficients[, -1])
      model$coefficients[, 1] <- 0
      
      initialvalues <- c(
        if ((npTM + sum(npEM) + npIPAll) > 0) {
          log(c(
            if (npTM > 0) model$transition_probs[paramTM],
            if (sum(npEM) > 0) {
              unlist(sapply(
                1:model$n_channels,
                function(x) model$emission_probs[[x]][paramEM[[x]]]
              ))
            },
            if (npIPAll > 0) {
              unlist(sapply(1:original_model$n_clusters, \(m) {
                if (npIP[m] > 0) original_model$initial_probs[[m]][paramIP[[m]]]
              }))
            }
          ))
        },
        model$coefficients[, -1]
      )
    } else {
      maxIP <- which.max(model$initial_probs)
      maxIPvalue <- model$initial_probs[maxIP]
      if (is.null(fixed_inits)) {
        initNZ <- model$initial_probs > 0
      } else {
        initNZ <- model$initial_probs > 0 & !fixed_inits
      }
      initNZ[maxIP] <- 0
      paramIP <- which(initNZ == 1)
      
      npIP <- length(paramIP)
      
      initialvalues <- c(if ((npTM + sum(npEM) + npIP) > 0) {
        log(c(
          if (npTM > 0) model$transition_probs[paramTM],
          if (sum(npEM) > 0) {
            unlist(sapply(
              1:model$n_channels,
              function(x) model$emission_probs[[x]][paramEM[[x]]]
            ))
          },
          if (npIP > 0) model$initial_probs[paramIP]
        ))
      })
    }
    
    
    scalerows <- function(x, fixed) {
      x[!fixed] <- x[!fixed] * (1 - sum(x[fixed])) / sum(x[!fixed])
      x
    }
    
    objectivef_mhmm <- function(pars, model, estimate = TRUE) {
      if (any(!is.finite(exp(pars))) && estimate) {
        if (need_grad) {
          return(list(objective = Inf, gradient = rep(-Inf, length(pars))))
        } else {
          Inf
        }
      }
      if (npTM > 0) {
        model$transition_probs[maxTM] <- maxTMvalue
        model$transition_probs[paramTM] <- exp(pars[1:npTM])
        zeros <- which(rowSums(model$transition_probs) == 0)
        diag(model$transition_probs)[zeros] <- 1
        if (is.null(fixed_transitions)) {
          model$transition_probs[] <-
            model$transition_probs / rowSums(model$transition_probs)
        } else {
          for (i in 1:nrow(model$transition_probs)) {
            cols <- fixed_transitions[i, ]
            if (any(cols)) {
              model$transition_probs[i, ] <- scalerows(model$transition_probs[i, ], cols)
            } else {
              model$transition_probs[i, ] <-
                model$transition_probs[i, ] / sum(model$transition_probs[i, ])
            }
          }
        }
      }
      if (sum(npEM) > 0) {
        if (is.null(constraints)) {
          for (i in 1:model$n_channels) {
            emissionArray[, 1:model$n_symbols[i], i][maxEM[[i]]] <- maxEMvalue[[i]]
            emissionArray[, 1:model$n_symbols[i], i][paramEM[[i]]] <-
              exp(pars[(npTM + 1 + c(0, cumsum(npEM))[i]):(npTM + cumsum(npEM)[i])])
            
            if (is.null(fixed_emissions)) {
              emissionArray[, 1:model$n_symbols[i], i] <-
                emissionArray[, 1:model$n_symbols[i], i] / rowSums(emissionArray[, 1:model$n_symbols[i], i, drop = FALSE])
            } else {
              for (j in 1:nrow(model$transition_probs)) {
                cols <- fixed_emissions[[i]][j, ]
                if (any(cols)) {
                  emissionArray[j, 1:model$n_symbols[i], i] <- scalerows(emissionArray[j, 1:model$n_symbols[i], i], cols)
                } else {
                  emissionArray[j, 1:model$n_symbols[i], i] <-
                    emissionArray[j, 1:model$n_symbols[i], i] / sum(emissionArray[j, 1:model$n_symbols[i], i, drop = FALSE])
                }
              }
            }
            
            emissionArray[, 1:model$n_symbols[i], i] <-
              emissionArray[, 1:model$n_symbols[i], i] / rowSums(emissionArray[, 1:model$n_symbols[i], i])
          }
        } else {
          emArray <- emissionArray[uniqs, , , drop = FALSE]
          for (i in 1:model$n_channels) {
            emArray[, 1:model$n_symbols[i], i][maxEM[[i]]] <- maxEMvalue[[i]]
            emArray[, 1:model$n_symbols[i], i][paramEM[[i]]] <-
              exp(pars[(npTM + 1 + c(0, cumsum(npEM))[i]):(npTM + cumsum(npEM)[i])])
            
            if (is.null(fixed_emissions)) {
              emArray[, 1:model$n_symbols[i], i] <-
                emArray[, 1:model$n_symbols[i], i] / rowSums(emArray[, 1:model$n_symbols[i], i, drop = FALSE])
            } else {
              for (j in 1:nrow(model$transition_probs)) {
                cols <- fixed_emissions[[i]][j, ]
                if (any(cols)) {
                  emArray[j, 1:model$n_symbols[i], i] <- scalerows(emArray[j, 1:model$n_symbols[i], i], cols)
                } else {
                  emArray[j, 1:model$n_symbols[i], i] <-
                    emArray[j, 1:model$n_symbols[i], i] / sum(emArray[j, 1:model$n_symbols[i], i, drop = FALSE])
                }
              }
            }
          }
          for (i in 1:n_states) {
            emissionArray[constraints == i, , ] <- rep(emArray[i, , ], each = sum(constraints == i))
          }
        }
      }
      for (m in 1:original_model$n_clusters) {
        if (npIP[m] > 0) {
          original_model$initial_probs[[m]][maxIP[[m]]] <- maxIPvalue[[m]]
          original_model$initial_probs[[m]][paramIP[[m]]] <- exp(pars[npTM + sum(npEM) + c(0, cumsum(npIP))[m] +
                                                                        1:npIP[m]])
          
          if (is.null(fixed_inits)) {
            original_model$initial_probs[[m]][] <- original_model$initial_probs[[m]] / sum(original_model$initial_probs[[m]])
          } else {
            original_model$initial_probs[[m]][] <- scalerows(original_model$initial_probs[[m]], fixed_inits[[m]])
          }
        }
      }
      model$initial_probs <- unlist(original_model$initial_probs)
      model$coefficients[, -1] <- pars[npTM + sum(npEM) + npIPAll + 1:npCoef]
      
      if (estimate) {
        if (!log_space) {
          if (need_grad) {
            objectivex(
              model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols,
              model$coefficients, model$X, model$n_states_in_clusters, threads
            )
          } else {
            -sum(logLikMixHMM(
              model$transition_probs, emissionArray, model$initial_probs, obsArray,
              model$coefficients, model$X, model$n_states_in_clusters, threads
            ))
          }
        } else {
          if (need_grad) {
            log_objectivex(
              model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols,
              model$coefficients, model$X, model$n_states_in_clusters, threads
            )
          } else {
            -sum(log_logLikMixHMM(
              model$transition_probs, emissionArray, model$initial_probs, obsArray,
              model$coefficients, model$X, model$n_states_in_clusters, threads
            ))
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
    
    objectivef_hmm <- function(pars, model, estimate = TRUE) {
      if (any(!is.finite(exp(pars))) && estimate) {
        if (need_grad) {
          return(list(objective = Inf, gradient = rep(-Inf, length(pars))))
        } else {
          Inf
        }
      }
      
      if (npTM > 0) {
        model$transition_probs[] <- original_trans
        model$transition_probs[paramTM] <- exp(pars[1:npTM])
        zeros <- which(rowSums(model$transition_probs) == 0)
        diag(model$transition_probs)[zeros] <- 1
        if (is.null(fixed_transitions)) {
          model$transition_probs[] <-
            model$transition_probs / rowSums(model$transition_probs)
        } else {
          for (i in 1:nrow(model$transition_probs)) {
            cols <- fixed_transitions[i, ]
            if (any(cols)) {
              model$transition_probs[i, ] <- scalerows(model$transition_probs[i, ], cols)
            } else {
              model$transition_probs[i, ] <-
                model$transition_probs[i, ] / sum(model$transition_probs[i, ])
            }
          }
        }
      }
      
      if (sum(npEM) > 0) {
        if (is.null(constraints)) {
          for (i in 1:model$n_channels) {
            emissionArray[, 1:model$n_symbols[i], i][maxEM[[i]]] <- maxEMvalue[[i]]
            emissionArray[, 1:model$n_symbols[i], i][paramEM[[i]]] <-
              exp(pars[(npTM + 1 + c(0, cumsum(npEM))[i]):(npTM + cumsum(npEM)[i])])
            
            if (is.null(fixed_emissions)) {
              emissionArray[, 1:model$n_symbols[i], i] <-
                emissionArray[, 1:model$n_symbols[i], i] / rowSums(emissionArray[, 1:model$n_symbols[i], i, drop = FALSE])
            } else {
              for (j in 1:nrow(model$transition_probs)) {
                cols <- fixed_emissions[[i]][j, ]
                if (any(cols)) {
                  emissionArray[j, 1:model$n_symbols[i], i] <- scalerows(emissionArray[j, 1:model$n_symbols[i], i], cols)
                } else {
                  emissionArray[j, 1:model$n_symbols[i], i] <-
                    emissionArray[j, 1:model$n_symbols[i], i] / sum(emissionArray[j, 1:model$n_symbols[i], i, drop = FALSE])
                }
              }
            }
          }
        } else {
          emArray <- emissionArray[uniqs, , , drop = FALSE]
          for (i in 1:model$n_channels) {
            emArray[, 1:model$n_symbols[i], i][maxEM[[i]]] <- maxEMvalue[[i]]
            emArray[, 1:model$n_symbols[i], i][paramEM[[i]]] <-
              exp(pars[(npTM + 1 + c(0, cumsum(npEM))[i]):(npTM + cumsum(npEM)[i])])
            
            if (is.null(fixed_emissions)) {
              emArray[, 1:model$n_symbols[i], i] <-
                emArray[, 1:model$n_symbols[i], i] / rowSums(emArray[, 1:model$n_symbols[i], i, drop = FALSE])
            } else {
              for (j in 1:nrow(model$transition_probs)) {
                cols <- fixed_emissions[[i]][j, ]
                if (any(cols)) {
                  emArray[j, 1:model$n_symbols[i], i] <- scalerows(emArray[j, 1:model$n_symbols[i], i], cols)
                } else {
                  emArray[j, 1:model$n_symbols[i], i] <-
                    emArray[j, 1:model$n_symbols[i], i] / sum(emArray[j, 1:model$n_symbols[i], i, drop = FALSE])
                }
              }
            }
          }
          for (i in 1:n_states) {
            emissionArray[constraints == i, , ] <- rep(emArray[i, , ], each = sum(constraints == i))
          }
        }
      }
      if (npIP > 0) {
        model$initial_probs[maxIP] <- maxIPvalue
        model$initial_probs[paramIP] <- exp(pars[(npTM + sum(npEM) + 1):(npTM + sum(npEM) + npIP)])
        if (is.null(fixed_inits)) {
          model$initial_probs[] <- model$initial_probs / sum(model$initial_probs)
        } else {
          model$initial_probs[] <- scalerows(model$initial_probs, fixed_inits)
        }
      }
      
      
      if (estimate) {
        if (!log_space) {
          if (need_grad) {
            objective(
              model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols, threads
            )
          } else {
            -sum(logLikHMM(
              model$transition_probs, emissionArray,
              model$initial_probs, obsArray, threads
            ))
          }
        } else {
          if (need_grad) {
            log_objective(
              model$transition_probs, emissionArray, model$initial_probs, obsArray,
              transNZ, emissNZ, initNZ, model$n_symbols, threads
            )
          } else {
            -sum(log_logLikHMM(
              model$transition_probs, emissionArray,
              model$initial_probs, obsArray, threads
            ))
          }
        }
      } else {
        if (sum(npEM) > 0) {
          if (!is.null(constraints)) {
            model$emission_probs <- original_emiss
          }
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
          lb <- c(
            rep(-25, length(initialvalues) - npCoef),
            rep(-150 / apply(abs(model$X), 2, max), model$n_clusters - 1)
          )
        } else {
          lb <- -25
        }
      }
      lb <- pmin(lb, 2 * initialvalues)
      if (missing(ub)) {
        if (mhmm) {
          ub <- c(
            rep(25, length(initialvalues) - npCoef),
            rep(150 / apply(abs(model$X), 2, max), model$n_clusters - 1)
          )
        } else {
          ub <- 25
        }
      }
      ub <- pmin(pmax(ub, 2 * initialvalues), 500)
      
      if (is.null(control_global$maxeval)) {
        control_global$maxeval <- 10000
      }
      if (is.null(control_global$maxtime)) {
        control_global$maxtime <- 60
      }
      if (is.null(control_global$algorithm)) {
        control_global$algorithm <- "NLOPT_GD_MLSL_LDS"
        if (is.null(control_global$local_opts)) {
          control_global$local_opts <- list(
            algorithm = "NLOPT_LD_LBFGS",
            ftol_rel = 1e-6, xtol_rel = 1e-4
          )
        }
      }
      need_grad <- grepl("NLOPT_GD_", control_global$algorithm)
      if (!is.null(constraints) && need_grad) {
        need_grad <- FALSE # don't use analytical gradients
        if (mhmm) {
          gr <- function(x, model, estimate) {
            nl.grad(x, objectivef_mhmm,
                    model = model, estimate = estimate
            )
          }
          globalres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_mhmm, lb = lb, ub = ub,
            eval_grad_f = gr,
            opts = control_global, model = model, estimate = TRUE, ...
          )
          model <- objectivef_mhmm(globalres$solution, model, FALSE)
        } else {
          gr <- function(x, model, estimate) {
            nl.grad(x, objectivef_hmm,
                    model = model, estimate = estimate
            )
          }
          globalres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_hmm, lb = lb, ub = ub,
            eval_grad_f = gr,
            opts = control_global, model = model, estimate = TRUE, ...
          )
          model <- objectivef_hmm(globalres$solution, model, FALSE)
        }
      } else {
        if (mhmm) {
          globalres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_mhmm, lb = lb, ub = ub,
            opts = control_global, model = model, estimate = TRUE, ...
          )
          model <- objectivef_mhmm(globalres$solution, model, FALSE)
        } else {
          globalres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_hmm, lb = lb, ub = ub,
            opts = control_global, model = model, estimate = TRUE, ...
          )
          model <- objectivef_hmm(globalres$solution, model, FALSE)
        }
      }
      
      if (globalres$status < 0) {
        warning_(paste("Global optimization terminated:", globalres$message))
      }
      initialvalues <- globalres$solution
      ll <- -globalres$objective
    } else {
      globalres <- NULL
    }
    
    if (local_step) {
      if (is.null(control_local$maxeval)) {
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
      
      need_grad <- grepl("NLOPT_LD_", control_local$algorithm)
      if ((!is.null(constraints) | fixed) && need_grad) {
        need_grad <- FALSE # don't use analytical gradients (do not take account constraints)
        if (mhmm) {
          gr <- function(x, model, estimate) {
            nl.grad(x, objectivef_mhmm,
                    model = model, estimate = estimate
            )
          }
          localres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_mhmm,
            eval_grad_f = gr,
            opts = control_local, model = model, estimate = TRUE, ...
          )
          model <- objectivef_mhmm(localres$solution, model, FALSE)
        } else {
          gr <- function(x, model, estimate) {
            nl.grad(x, objectivef_hmm,
                    model = model, estimate = estimate
            )
          }
          localres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_hmm,
            eval_grad_f = gr,
            opts = control_local, model = model, estimate = TRUE, ...
          )
          model <- objectivef_hmm(localres$solution, model, FALSE)
        }
      } else {
        if (mhmm) {
          localres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_mhmm,
            opts = control_local, model = model, estimate = TRUE, ...
          )
          model <- objectivef_mhmm(localres$solution, model, FALSE)
        } else {
          localres <- nloptr(
            x0 = initialvalues, eval_f = objectivef_hmm,
            opts = control_local, model = model, estimate = TRUE, ...
          )
          model <- objectivef_hmm(localres$solution, model, FALSE)
        }
      }
      if (localres$status < 0) {
        warning_(paste("Local optimization terminated:", localres$message))
      }
      ll <- -localres$objective
    } else {
      localres <- NULL
    }
    
    if (mhmm) {
      rownames(model$coefficients) <- colnames(model$X)
      colnames(model$coefficients) <- model$cluster_names
    }
  } else {
    globalres <- localres <- NULL
  }
  
  if (model$n_channels == 1) {
    model$emission_probs <- model$emission_probs[[1]]
  }
  
  if (mhmm) {
    model <- spread_models(model)
    for (i in 1:model$n_clusters) {
      dimnames(model$transition_probs[[i]]) <- dimnames(original_model$transition_probs[[i]])
      for (j in 1:model$n_channels) {
        dimnames(model$emission_probs[[i]][[j]]) <- dimnames(original_model$emission_probs[[i]][[j]])
      }
    }
  }
  
  suppressWarnings(try(model <- trim_model(model, verbose = FALSE), silent = TRUE))
  list(
    model = model,
    logLik = ll, em_results = resEM[c("logLik", "iterations", "change", "best_opt_restart")],
    global_results = globalres, local_results = localres, call = match.call()
  )
}
