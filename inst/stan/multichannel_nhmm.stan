functions {
  #include /include/array_additions.stan
  #include /include/loglik_multichannel.stan
}
data {
  #include /include/data_multichannel.stan
  #include /include/data_covariates.stan
}
transformed data {
  #include /include/transformed_data.stan
}
parameters {
  #include /include/parameters_multichannel.stan
}
transformed parameters {
  #include /include/priors_multichannel.stan
  real log_lik = loglik_mc(beta_i_raw, beta_s_raw, beta_o_raw, obs, M, T, 
  N_sample, ids, X_s, X_o, X_i);
}
model {
  target += prior;
  target += log_lik;
}
generated quantities {
  real ploglik_N = prior + log_lik;
  if (compute_ploglik_N == 1) {
    ploglik_N = prior + loglik_mc(beta_i_raw, beta_s_raw, beta_o_raw, obs, 
    M, T, N, linspaced_int_array(N, 1, N), X_s, X_o, X_i);
  }
}
