functions {
  #include /include/loglik.stan
  #include /include/array_additions.stan
}
data {
  #include /include/data_multichannel.stan
  #include /include/data_covariates.stan
}
transformed data {
  #include /include/transformed_data.stan
  matrix[M, K_o] zeros_mat_M_K_o = rep_matrix(0, M, K_o);
  array[C + 1] int cumsum_M = append_array({0}, cumulative_sum(n_M));
}
parameters {
  #include /include/parameters_multichannel.stan
}
transformed parameters {
  real log_lik;
  #include /include/transformed_parameters_multichannel.stan
  #include /include/priors_multichannel.stan
  {
    vector[N] ll;
    vector[S] log_Pi;
    array[C] matrix[S, M + 1] log_B;
    matrix[S, T] log_py;
    array[T] matrix[S, S] log_A;
    vector[S] zeros = rep_vector(0, S);
    for(i in 1:N) {
      #include /include/model_pi_varying.stan
      for(t in 1:T) {
        #include /include/model_A_varying.stan
        #include /include/model_B_varying_multichannel.stan
        log_py[, t] = zeros;
        for(c in 1:C) {
          log_py[, t] += log_B[c, , obs[c, t, i]];
        }
      }
      ll[i] = loglik(log_Pi, log_A, log_py);
    }
    log_lik = sum(ll);
  }
}
model {
  target += prior;
  target += log_lik;
}
