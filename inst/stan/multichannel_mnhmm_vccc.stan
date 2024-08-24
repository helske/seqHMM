functions {
  #include /include/loglik.stan
  #include /include/array_additions.stan
}
data {
  #include /include/data_multichannel.stan
  #include /include/data_covariates.stan
  int<lower=2> D;
  int<lower=1> K_d; // number of covariates for initial probabilities (including the intercept)
  array[N] vector[K_d] X_d; // covariates for initial probabilities (including the intercept)
}
transformed data {
  #include /include/transformed_data_multichannel.stan
  row_vector[K_d] zeros_row_K_d = rep_row_vector(0, K_d);
}
parameters {
  #include /include/parameters_multichannel_mixture.stan
}
transformed parameters {
  real log_lik;
  #include /include/transformed_parameters_multichannel_mixture.stan
  #include /include/priors_multichannel_mixture.stan
  {
    matrix[N, D] ll;
    vector[S] log_Pi;
    array[C] matrix[S, max_M + 1] log_B;
    matrix[S, T] log_py;
    matrix[S, S] log_A;
    for(d in 1:D) {
      #include /include/model_A_constant_mixture.stan
      #include /include/model_B_constant_multichannel_mixture.stan
      for(i in 1:N) {
        #include /include/model_pi_varying_mixture.stan
        for(t in 1:T) {
          log_py[, t] = zeros_S;
          for(c in 1:C) {
            log_py[, t] += log_B[c, , obs[c, t, i]];
          }
        }
        ll[i, d] = loglik(log_Pi, log_A, log_py);
      }
    }
    #include /include/model_theta_constant_mixture.stan
    vector[N] ll_i;
    for (i in 1:N) {
      ll_i[i] = log_sum_exp(log_omega + ll[i, ]');
    }
    log_lik = sum(ll_i);
  }
}
model {
  target += prior;
  target += log_lik;
}
