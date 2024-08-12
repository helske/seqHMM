int<lower=1> K_o; // number of covariates for emission probabilities
array[T, N] vector[K_o] X_o; // covariates for emissions (including the intercept)
