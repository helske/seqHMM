//log-sum-exp trick
#include "seqHMM.h"

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
  unsigned int maxi = 0;
  LDOUBLE maxv = x.max(maxi);
  if (!(maxv > -arma::math::inf())) {
    return -arma::math::inf();
  }
  LDOUBLE cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -arma::math::inf())) {
      cumsum += EXPL(x(i) - maxv);
    }
  }
  
  return maxv + log1p(cumsum);
}
