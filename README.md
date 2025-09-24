[![R-CMD-check](https://github.com/helske/seqHMM/workflows/R-CMD-check/badge.svg)](https://github.com/helske/seqHMM/actions)
[![Codecov test coverage](https://codecov.io/gh/helske/seqHMM/branch/main/graph/badge.svg)](https://app.codecov.io/gh/helske/seqHMM?branch=main)
[![cran version](https://www.r-pkg.org/badges/version/seqHMM)](https://CRAN.R-project.org/package=seqHMM)
[![downloads](https://cranlogs.r-pkg.org/badges/seqHMM)](https://cranlogs.r-pkg.org/badges/seqHMM)

seqHMM: Hidden Markov Models for Life Sequences and Other Multivariate, Multichannel Categorical Time Series
====================================================================================================

The `seqHMM` package is designed for fitting hidden (latent) Markov models (HMMs) and their variations for social sequence data and other categorical sequence data (e.g. categorical time series and panel data). Restricted and extended variants include mixture HMMs, Markov models and their mixtures,  latent class models, non-homogeneous hidden Markov models (NHMMs) and their mixtures, and feedback-augmented hidden Markov models (FAN-HMMs). 

The package supports models for one or multiple subjects with one or multiple parallel outcome sequences (channels). External covariates can be added to explain cluster membership in mixture models, and NHMMs and their variants support covariates in initial, transition and emission matrices as well.

Maximum likelihood estimation via EM algorithm and direct numerical maximization with analytical gradients is supported. All main algorithms are written in C++. Parallel computation is available via OpenMP (pre-2.0.0 version models) and [`future`](https://future.futureverse.org/) (via parallel multistart estimation with random initial values).

When using the package in publications, please cite:

Helske, Satu and Helske, Jouni (2019). Mixture hidden Markov models for sequence data: the seqHMM package in R. *Journal of Statistical Software, 88*(3). [doi:10.18637/jss.v088.i03](https://dx.doi.org/10.18637/jss.v088.i03).

Helske, Jouni (2025). Feedback-augmented Non-homogeneous Hidden Markov Models for Longitudinal Causal Inference. *arXiv preprint*. [doi:10.48550/arXiv.2503.16014](https://doi.org/10.48550/arXiv.2503.16014).

If you find bugs, please add a new issue here in GitHub. You can also contact Satu Helske (firstname.lastname@utu.fi) or Jouni Helske (firstname.lastname@iki.fi). We would be happy to hear your feedback.

The package is available on CRAN. Install it via

```R
install.packages("seqHMM")
```
If you want to try the development version of the `seqHMM` package, install it from the [R-universe](https://helske.r-universe.dev/seqHMM):

```R
install.packages("seqHMM", repos = "https://helske.r-universe.dev")
```

