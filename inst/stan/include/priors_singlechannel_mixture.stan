// priors for (very weak) regularisation
real prior = 0;
if (penalize == 1) {
  for (d in 1:D){
    prior += normal_lpdf(to_vector(beta_i_raw[d]) | 0, penalty);
    for(s in 1:S) {
      prior += normal_lpdf(to_vector(beta_s_raw[d, s]) | 0, penalty);
      prior += normal_lpdf(to_vector(beta_o_raw[d, s]) | 0, penalty);
    }
  }
  prior += normal_lpdf(to_vector(theta_raw) | 0, penalty);
}
