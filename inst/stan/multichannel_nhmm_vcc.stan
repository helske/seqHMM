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
    array[C] matrix[S, M] log_B;
    matrix[S, T] log_py;
    matrix[S, S] log_A;
    #include /include/model_A_constant.stan
    #include /include/model_B_constant_multichannel.stan
    vector[S] zeros = rep_vector(0, S);
    for(i in 1:N) {
      #include /include/model_pi_varying.stan
      for(t in 1:T) {
        log_py[, t] = zeros;
        for(c in 1:C) {
          if(obs[c, t, i] > 0) {
            log_py[, t] += log_B[c, , obs[c, t, i]];
          }
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
