  // priors for (very weak) regularisation
  real prior = normal_lpdf(to_vector(beta_i_raw) | 0, 5);
  for(s in 1:S) {
    prior = normal_lpdf(to_vector(beta_s_raw[s]) | 0, 5);
  }
  prior = normal_lpdf(beta_o_raw | 0, 5);
