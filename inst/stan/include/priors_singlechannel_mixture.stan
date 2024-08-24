  // priors for (very weak) regularisation
  real prior = 0;
  for (d in 1:D){
    prior += normal_lpdf(to_vector(beta_i_raw[d]) | 0, 5);
    for(s in 1:S) {
      prior += normal_lpdf(to_vector(beta_s_raw[d, s]) | 0, 5);
      prior += normal_lpdf(to_vector(beta_o_raw[d, s]) | 0, 5);
    }
  }
