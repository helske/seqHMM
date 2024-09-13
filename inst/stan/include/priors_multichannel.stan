real prior = 0;
if (penalize == 1) {
  // priors for (very weak) regularisation
  prior += normal_lpdf(to_vector(beta_i_raw) | 0, penalty);
  for(s in 1:S) {
    prior += normal_lpdf(to_vector(beta_s_raw[s]) | 0, penalty);
  }
  prior += normal_lpdf(beta_o_raw | 0, penalty);
}
