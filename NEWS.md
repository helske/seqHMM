seqHMM 2.1.1
==============
  * Fixed a bug in `coef.mnhmm` which resulted an error on extracting 
  confidence intervals for emission probabilities.

seqHMM 2.1.0
==============

  * Identifiability check via QR decomposition for the design matrices can now 
  be controlled via argument `check_rank`. Default is `NULL`, which will 
  evaluate as `TRUE` if the number of sequences is at most 1000. Also, even if 
  checks are done, they now use precomputed QR instead of doing it again.
  * Added dependency on package `collapse`, which resulted in much more 
  efficient construction of non-homogeneous HMM objects.
  * Removed `data.table` progress printing from internal functions. Because of 
  this recent feature, `data.table` version requirement is increased.
  * Fixed an issue with dropping the first time point in case where this lead 
  to an unused factor level and subsequent error about multicollinearity.
  * Fixed the error messages regarding the missing values in covariates.
  * Fixed a typo in `simulate_nhmm` which caused an error in accessing 
  elements of initial values.
  * Fixed the bug in `simulate_nhmm` and `simulate_mnhmm` which caused the 
  covariate matrices of the returned model to be incorrect. This had no effect 
  on the simulated observations and states though.
  * Unused factor levels are now not dropped in the simulation functions, 
  whereas whether to drop them in estimation functions can be controlled via 
  argument `drop_levels` (which is `TRUE` by default).
  * Fixed the column name from `state_to` to `state` in the marginal state 
  predictions output of `predict()`.
  * Fixed the check for valid `ids` argument in `stacked_sequence_plot`.
  * Fixed the usage of the `group` argument in `stacked_sequence_plot`.
  * Reverted the default argument `log_space` of `fit_model` 
    and other relevant functions back to `FALSE`, as was the case before 
    package version 2.0.0. Algorithms of NHMM models also now use scaling 
    instead of log-space computations for improved speed.
  * Improved numerical stability and computational efficiency of several 
    C++ functions related to non-homogenous models. 
  * Resaved the `leaves` data as normal `data.frame` (instead of `tibble`).
  
seqHMM 2.0.0
==============
  * Added support for non-homogeneous HMMs (NHMMs) where initial, transition, 
  and emission probabilities can depend on individual-specific covariates. For 
  mixture NHMMs cluster probabilities can also depend on covariates.
  * Added support for feedback-augmented NHMMs where responses and states can 
  depend also on the past responses.
  * Added `bootstrap_coef` bootstrapping coefficients of NHMMs.
  * Added a function `predict`, which can be used to compute average marginal 
  predictions for NHMMs, which can be interpreted as average causal effects 
  under suitable assumptions (https://arxiv.org/abs/2503.16014).
  * Added `get_*` functions for obtaining initial, transition, emission and 
  cluster probabilities and their (conditional) marginals.
  * Rewrote all post-estimation functions such as `hidden_paths` and 
  `forward_backward` so that they return `data.table`.
  * Added functions `stslist_to_data` and `data_to_stslist` to convert between 
  data frames and TraMineR's `stslist` objects (created by `seqdef`).
  * Rewrote sequence visualization functions using `ggseqplot` and `patchwork` 
  packages. Old plotting functions are deprecated and will be removed in the 
  future.
  * Warning and error messages were rewritten using `cli` package.
  * Removed the function `estimate_coef`. A portion of the code was accidentally 
  commented out, rendering the function non-functional for several years. 
  Rather than correcting the code, the function was removed as it was deemed 
  unnecessary.
  * Added automatic tests using `testthat` package.
  * Internally switched from Rd syntax to Markdown in Roxygen documentation.

seqHMM 1.2.6 (Release date: 2023-06-07)
==============
  * Armadillo solver now fails directly instead of trying approximate solution 
    in case of numerical issues when estimating the regression coefficients.
  * `build_mm` now correctly runs EM algorithm only in case of real missing 
    observations in data, i.e. it correctly handles void values.
    
seqHMM 1.2.5 (Release date: 2023-06-13)
==============
  * New methods `state_names` and `cluster_names` for getting and setting state 
    and cluster names of the model objects.
  * Fixed a bug in `mssplot` which caused error with `respect_void = TRUE`.
  * Fixed a bug in `simulate_mhmm` which resulted NA colnames in the simulated 
    `stslist` objects.
  * Modified internals of `ssp` to be compatible with the latest versions of 
    TraMineR.
  * The `data` argument is no longer optional if `formula` is used in 
    `build_mhmm`, `build_mmm` and `build_lcm` functions.
  * Added few unit tests using testthat.

seqHMM 1.2.4 (Release date: 2023-01-09)
==============
 * Fixed the issue with the initial transition matrix construction in build_mm 
   when there is a symbol which is only present in last time points.
 * Related to above, if a state x with no transitions from is encountered during 
   the EM algorithm, the corresponding row of the transition matrix A is then 
   normalized so that is such a A[x,x]=1 i.e. the state is set as absorbing.
 
seqHMM 1.2.3 (Release date: 2022-12-13)
==============
 * Changed the internal seqplot axes argument to xaxis as axes is now 
   deprecated in TraMineR.
   
seqHMM 1.2.2 (Release date: 2022-12-01)
==============
* Fixed the handling of models with more than 200 symbols in mc_to_sc, 
  mc_to_sc_data, and plotting functions.
* Removed duplicate xaxis argument in interal SSPlotter function.
  
seqHMM 1.2.1-1 (Release date: 2022-05-24)
==============
* Added argument respect_void to hidden_paths function leading to the 
  propagation of void values of observed sequences to the hidden state sequences.
* Fixed the scaling of initial probabilities to 1 in build_mm.
* Fixed the computation of the initial state probability vector for Markov 
  models in case of missing data.

seqHMM 1.2.0 (Release date: 2021-10-18)
==============
* seqHMM now supports fixed values in initial, transition, and emission 
  probabilities when performing model estimation.
  

seqHMM 1.1.1 (Release date: 2021-08-13)
==============
* Fixed a a case with the equality constraints which resulted an error during 
  the estimation in case of (numerically) nonfinite parameter value.
* Function fit_model now checks whether the initial model produces finite 
  likelihood before proceeding further.
  
seqHMM 1.1.0 (Release date: 2021-06-18)
==============
* Added a feature which allows equality constraints for emission probability distributions.
* Fixed the documentation of biofam3c data.
* Added a check for degenerate model for 'fit_model' 
  where one channels contains only missing values.

seqHMM 1.0.14 (Release date: 2019-10-21)
==============
* Added an option for defining legend labels and colours manually in plot.hmm.

seqHMM 1.0.13 (Release date: 2019-06-11)
==============
* Fixed a bug when calling plot.mhmm inside a function. Thanks to Ellen Graham for 
  catching the issue.
  
seqHMM 1.0.12 (Release date: 2019-04-11)
==============
* Fixed an OpenMP sharing issue due to upcoming GCC9.

seqHMM 1.0.11 (Release date: 2019-04-09)
==============
* All functions relying on C++ are now slightly faster due to 
  more efficient data transformations from factors to integers.
* Fixed an OpenMP sharing issue due to upcoming GCC9.

seqHMM 1.0.10 (Release date: 2019-01-25)
==============
* Fixed iteration counter for EM algorithm, now maxeval=1 actually does 
  one iteration and not two. Similarly fixed the the corresponding printing of 
  intermediate results.
* build_hmm and build_hmm now ignore n_states argument if any of the initial values 
  are provided.
* Updated some formulas in the algorithms vignette to correspond to current forward-backward   variant.
* Updated citation info due to the publication in JSS.

seqHMM 1.0.9 (Release date: 2018-11-06)
==============
* Fixed a bug in backward and posterior probabilities of multivariate models stemming 
  from changes in version 1.0.7, causing the incorrect scaling factors. Note that issue
  was only in forward_backward and posterior_probs functions with log_scale = FALSE, 
  internally the package still used correct implementations when calling fit_model.
  Thanks for Silvia Bacci for noticing the issue.
* Fixed OpenMP flags in Makevars.

seqHMM 1.0.8-1 (Release date: 2018-05-03)
==============
* Fixed encoding issues in references.
* Updated author affiliations.

seqHMM 1.0.8 (Release date: 2017-11-07)
==============

* Fixed a bug in forward_backward function, which did not provide correct 
  scaling factors for last time point and thus the backward variables weren't F
  scaled correctly.
* Related to above bug, posterior_probs now provides proper probabilities 
  between 0 and 1.
* Fixed a bug in ssp; caused an error for sortv = "hidden.paths" when hidden 
  paths were not provided even though x was an hmm object.
* Changed argument withlegend to with.legend due to changes in the TraMineR package. 
* Related to the TraMineR update, fixed warnings given by ssp functions.
* Fixing model plots after an update in the igraph package.
  
seqHMM 1.0.7 (Release date: 2017-04-04)
==============

* Added supplementary vignettes for visualization, estimation, and theoretical 
  background.
* Corrected a "feature" of the scaled forward-backward algorithm which caused 
  potential numerical issues in backward probabilities. More specifically, 
  previously the scaled backward variable at time t was scaled with the scaling 
  factors c_[t+1] of the forward variable, instead of c_t. In typical 
  applications this did not cause any problems, but some models which previously 
  caused problems can work now without resorting to the log-space algorithm.
* Scaling factors used in forward-backward algorithms are now stored as 1/c_t,
  where c_t are the scaling factors from older versions of seqHMM.
* The build_mm function now automatically estimates model parameters from the 
  observed initial state distribution and transition counts.
* Added automatic starting values for build_hmm, build_mhmm, build_mmm, and 
  build_lcm.
  
Bug fixes:

* Fixed a bug from version 1.0.6 where the EM algorithm with restarts did not 
  work as efficiently as intended due to changes in initializing emission matrices. 
  More specifically, the previous version used Rcpp sugar which does not use deep 
  copies in case of arma::cubes, which resulted the subsequent EM runs to depend 
  on the first optimum (affecting especially the identification of non-zero 
  emission probabilities).
* Fixed a bug in the ssp function where tlim was taken account of after sorting 
  sequences and computing hidden paths (now cases are chosen before other actions).


seqHMM 1.0.6 (Release date: 2016-08-01)
==============

* Argument diag_c in simulate_transition_probs is now used also in cases where
  left_right = FALSE.
* Adjusted reltol and maxeval values for EM algorithm. Now reltol is 1e-10
  (previously 1e-12), and the reltol and maxeval values for restarts are by
  default taken from the initial EM algorithm (previously reltol was 1e-8 and
  maxeval = 100 for restarts).
* Fixed hidden states labels for ssp functions (previously always used the 
  default values).

Bug fixes:

* The ssplot function assigned wrong colors for hidden states in cases where
  the state names were not alphabetically ordered. The performance of the
  function was also improved by removing extra calls to seqdef.
* Changing the missing.color argument did not work in legends of ssp, ssplot, and
  mssplot.
* The mssplot function now works with unique hidden state names (problem
  occurred e.g. with latent class models).
* The mssplot function with sortv = "mds.hidden" produced strange errors when
  plotting clusters with one hidden state. Now automatically uses "mds.obs" in
  such cases.

seqHMM 1.0.5 (Release date: 2016-02-24)
==============

Bug fixes:

* The mssplot function now uses hidden paths instead of posterior probabilities to determine the most probable cluster for each subject (previous solution caused errors when posterior probabilities suggested a different cluster than hidden paths).
* In mssplot, removed a misplaced tlim which caused a warning when plotting state distributions of hidden paths.
* The gridplot function now uses with.missing.legend also with combined legends.



seqHMM 1.0.4 (Release date: 2016-01-14)
==============

* Added examples for build_mmm.
* Added more space for main titles in plot.mhmm.
* Improved documentation.
* Added tlim in ssp functions.

Bug Fixes:

* Corrected mc_to_sc for single state models due to dimension dropping.
* Corrected a bug in plot.hmm when layout = "vertical".
* Fixed legend layout in plot.hmm.
* Wrong nobs and df attributes in mc_to_sc.
* Wrong number of sequences to ssp titles when tlim is used.
* EM with HMM using log-space produced error due to missing element in output of EM.


seqHMM 1.0.3-1 (Release date: 2015-12-29)
==============

As requested by CRAN, changed donttest example of interactive plotting to conditional block depending
on whether session is interactive or not.

seqHMM 1.0.3 (Release date: 2015-12-23)
==============

Corrected dependency on R 3.2.0 due to lengths function.

Bug Fixes:

* Corrected a bug which caused fit_model to stop if restarted EM failed.
* vcov.mhmm produced errors in valgrind, corrected issue by replacing Armadillo's inv_sympd function with inv.
* Corrected a bug relating to colorpalette in mc_to_sc function.

Performance improvements:

* Slight performance improvement in all functions by tweaking the usage of armadillo constructors.

seqHMM 1.0.2-1 (Release date: 2015-12-19)
==============

First version on CRAN.
