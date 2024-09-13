functions {
  #include /include/array_additions.stan
  #include /include/loglik_singlechannel_mixture.stan
}
data {
  #include /include/data_singlechannel.stan
  #include /include/data_covariates.stan
  int<lower=2> D;
  int<lower=1> K_d; // number of covariates for initial probabilities (including the intercept)
  array[N] vector[K_d] X_d; // covariates for initial probabilities (including the intercept)
}
transformed data {
  #include /include/transformed_data.stan
}
parameters {
  #include /include/parameters_singlechannel_mixture.stan
}
transformed parameters {
  #include /include/priors_singlechannel_mixture.stan
  real log_lik = loglik_sc_mix(beta_i_raw, beta_s_raw, beta_o_raw, theta_raw, 
  obs, M, T, N_sample, ids, X_s, X_o, X_i, X_d);
}
model {
  target += prior;
  target += log_lik;
}
generated quantities {
  real ploglik_N = prior + log_lik;
  if (compute_ploglik_N == 1) {
    ploglik_N = prior + loglik_sc_mix(beta_i_raw, beta_s_raw, 
    beta_o_raw, theta_raw, obs, M, T, N, linspaced_int_array(N, 1, N), X_s, X_o, 
    X_i, X_d);
  }
}
