functions {
  #include /include/array_additions.stan
  #include /include/loglik_singlechannel.stan
}
data {
  #include /include/data_singlechannel.stan
  #include /include/data_covariates.stan
}
parameters {
  #include /include/parameters_singlechannel.stan
}
transformed parameters {
  #include /include/priors_singlechannel.stan
  real log_lik = loglik_sc(beta_i_raw, beta_s_raw, beta_o_raw, obs, M, T, 
  N_sample, ids, X_s, X_o, X_i);
}
model {
  target += prior;
  target += log_lik;
}
generated quantities {
  real ploglik_N = prior + loglik_sc(beta_i_raw, beta_s_raw, beta_o_raw, obs, 
  M, T, N, linspaced_int_array(N, 1, N), X_s, X_o, X_i);
}
