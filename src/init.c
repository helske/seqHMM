#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP seqHMM_EM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_EMx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_estimate_coefs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_forwardbackward(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_forwardbackwardx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_EM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_EMx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_forwardbackward(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_forwardbackwardx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_logLikHMM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_logLikMixHMM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_logLikHMM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_logLikMixHMM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_objective(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_log_objectivex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_logSumExp(SEXP);
extern SEXP seqHMM_objective(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_objectivex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_varcoef(SEXP, SEXP);
extern SEXP seqHMM_viterbi(SEXP, SEXP, SEXP, SEXP);
extern SEXP seqHMM_viterbix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"seqHMM_EM",                   (DL_FUNC) &seqHMM_EM,                    9},
  {"seqHMM_EMx",                  (DL_FUNC) &seqHMM_EMx,                  12},
  {"seqHMM_estimate_coefs",       (DL_FUNC) &seqHMM_estimate_coefs,       12},
  {"seqHMM_forwardbackward",      (DL_FUNC) &seqHMM_forwardbackward,       6},
  {"seqHMM_forwardbackwardx",     (DL_FUNC) &seqHMM_forwardbackwardx,      9},
  {"seqHMM_log_EM",               (DL_FUNC) &seqHMM_log_EM,                9},
  {"seqHMM_log_EMx",              (DL_FUNC) &seqHMM_log_EMx,              12},
  {"seqHMM_log_forwardbackward",  (DL_FUNC) &seqHMM_log_forwardbackward,   6},
  {"seqHMM_log_forwardbackwardx", (DL_FUNC) &seqHMM_log_forwardbackwardx,  9},
  {"seqHMM_logLikHMM",            (DL_FUNC) &seqHMM_logLikHMM,             5},
  {"seqHMM_logLikMixHMM",         (DL_FUNC) &seqHMM_logLikMixHMM,          8},
  {"seqHMM_log_logLikHMM",        (DL_FUNC) &seqHMM_log_logLikHMM,         5},
  {"seqHMM_log_logLikMixHMM",     (DL_FUNC) &seqHMM_log_logLikMixHMM,      8},
  {"seqHMM_log_objective",        (DL_FUNC) &seqHMM_log_objective,         9},
  {"seqHMM_log_objectivex",       (DL_FUNC) &seqHMM_log_objectivex,       12},
  {"seqHMM_logSumExp",            (DL_FUNC) &seqHMM_logSumExp,             1},
  {"seqHMM_objective",            (DL_FUNC) &seqHMM_objective,             9},
  {"seqHMM_objectivex",           (DL_FUNC) &seqHMM_objectivex,           12},
  {"seqHMM_varcoef",              (DL_FUNC) &seqHMM_varcoef,               2},
  {"seqHMM_viterbi",              (DL_FUNC) &seqHMM_viterbi,               4},
  {"seqHMM_viterbix",             (DL_FUNC) &seqHMM_viterbix,              7},
  {NULL, NULL, 0}
};

void R_init_seqHMM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}