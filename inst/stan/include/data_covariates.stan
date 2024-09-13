  int<lower=1> K_i; // number of covariates for initial probabilities (including the intercept)
  array[N] vector[K_i] X_i; // covariates for initial probabilities (including the intercept)
  int<lower=1> K_s; // number of covariates for transition probabilities (including the intercept)
  array[max_T, N] vector[K_s] X_s; // covariates for transitions (including the intercept)
  int<lower=1> K_o; // number of covariates for emission probabilities (including the intercept)
  array[max_T, N] vector[K_o] X_o; // covariates for emissions (including the intercept)
  int<lower=0,upper=1> penalize; // penalize the likelihood of the model
  real<lower=0> penalty; // standard deviation of the priors
