#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _cPCG_cgsolve(SEXP, SEXP, SEXP, SEXP);
extern SEXP _cPCG_icc(SEXP);
extern SEXP _cPCG_pcgsolve(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_cPCG_cgsolve",  (DL_FUNC) &_cPCG_cgsolve,  4},
    {"_cPCG_icc",      (DL_FUNC) &_cPCG_icc,      1},
    {"_cPCG_pcgsolve", (DL_FUNC) &_cPCG_pcgsolve, 5},
    {NULL, NULL, 0}
};

void R_init_cPCG(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}