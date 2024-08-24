        for(s in 1:S) {
          log_A[t, s, ] = log_softmax(beta_s[d, s] * X_s[t, i])';
        }
