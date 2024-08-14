//log-sum-exp trick and softmax functions
#include "logsumexp.h"

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  LDOUBLE maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  LDOUBLE cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += EXPL(x(i) - maxv);
    }
  }
  
  return maxv + log1p(cumsum);
}

// [[Rcpp::export]]
arma::vec softmax(const arma::vec& x) {
  double x_max = arma::max(x);
  arma::vec result = arma::exp(x - x_max);
  result = result / sum(result);
  return result;
}
