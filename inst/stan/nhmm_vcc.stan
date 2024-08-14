functions {
  #include /include/loglik.stan
  #include /include/array_additions.stan
}
data {
  #include /include/data_singlechannel.stan
  #include /include/data_covariates.stan
}
transformed data {
  #include /include/transformed_data.stan
}
parameters {
  #include /include/parameters_singlechannel.stan
}
transformed parameters {
  real log_lik;
  #include /include/transformed_parameters_singlechannel.stan
  #include /include/priors_singlechannel.stan
  {
    vector[N] ll;
    vector[S] log_Pi;
    matrix[S, M] log_B;
    matrix[S, T] log_py;
    matrix[S, S] log_A;
    #include /include/model_A_constant.stan
    #include /include/model_B_constant.stan
    for(i in 1:N) {
      #include /include/model_pi_varying.stan
      for(t in 1:T) {
        if(obs[t, i]  == 0) {
          log_py[, t] = rep_vector(0, S);
        } else {
          log_py[, t] = log_B[, obs[t, i]];
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
