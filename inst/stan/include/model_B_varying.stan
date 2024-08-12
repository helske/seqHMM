for(s in 1:S) {
  log_B[s, ] = log_softmax(beta_o[s] * X_o[t, i])';
}
