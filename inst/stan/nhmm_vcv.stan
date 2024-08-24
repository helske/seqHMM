functions {
  #include /include/loglik.stan
  #include /include/array_additions.stan
}
data {
  #include /include/data_singlechannel.stan
  #include /include/data_covariates.stan
}
transformed data {
  #include /include/transformed_data_singlechannel.stan
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
    matrix[S, M + 1] log_B;
    matrix[S, T] log_py;
    matrix[S, S] log_A;
    #include /include/model_A_constant.stan
    for(i in 1:N) {
     log_Pi = log_softmax(beta_i * X_i[i]);
      for(t in 1:T) {
        #include /include/model_B_varying.stan
        log_py[, t] = log_B[, obs[t, i]];
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
