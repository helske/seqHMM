for (c in 1:C) {
  for (s in 1:S) {
    log_B[c, s, ] = log(B[c, s]');
  }
}
