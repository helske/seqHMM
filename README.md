seqHMM: Hidden Markov Models for Life Sequences and Other Multivariate, Multichannel Categorical Time Series
====================================================================================================

Package seqHMM is designed for inference of hidden Markov models
    where both the hidden state space and the symbol space of observations is
    discrete, and observations consists of multiple sequences possibly with
    multiply channels such as life calendar data with multiple life domains.
    Maximum likelihood estimation via direct numerical maximization with 
    analytical gradients and EM algorithm is supported. All main algorithms 
    are written in C++.

Package is still under heavy development (see details below), and should be available at CRAN in early 2015.

If you want to try it out, you can install it via devtools package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/seqHMM")
```

What is still missing:

<ul>
 <li>Pretty plotting functions via TraMiner are under active development</li>
 <li>Function for computing posterior probabilities is currently missing (not ported from Fortran version yet, just laziness)</li>
 <li>Support for covariates for initial states (relates to _working_ paper by Helske&Helske et al. dealing with HMM clustering)</li>
</ul> 



