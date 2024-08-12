functions {
  #include /include/loglik.stan
}
data {
  #include /include/data_multichannel.stan
}
parameters {
  #include /include/parameters_pi_constant.stan
  #include /include/parameters_A_constant.stan
  #include /include/parameters_B_constant_multichannel.stan
}
transformed parameters {
  real log_lik_;
  real prior_ = 0.0;
  #include /include/prior_pi_constant.stan
  #include /include/prior_A_constant.stan
  #include /include/prior_B_constant_multichannel.stan
  {
    vector[N] ll;
    vector[S] log_Pi;
    array[C] matrix[S, M] log_B;
    matrix[S, T] log_py;
    matrix[S, S] log_A;
    #include /include/model_pi_constant.stan
    #include /include/model_A_constant.stan
    #include /include/model_B_constant_multichannel.stan
    vector[S] zeros = rep_vector(0, S);
    for(i in 1:N) {
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
    log_lik_ = sum(ll);
  }
}
model {
  target += prior_;
  target += log_lik_;
}
generated quantities {
  real prior = prior_;
  real log_lik = log_lik_;
}
