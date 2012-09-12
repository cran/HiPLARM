#ifndef PLASMA_POMATRIX_H
#define PLASMA_POMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP plasma_dpoMatrix_rcond(SEXP obj, SEXP type);
SEXP plasma_dpoMatrix_solve(SEXP a);
SEXP plasma_dpoMatrix_matrix_solve(SEXP a, SEXP b);
SEXP plasma_dpoMatrix_dgeMatrix_solve(SEXP a, SEXP b);
SEXP plasma_dpoMatrix_chol(SEXP x);
double get_norm_sy(SEXP obj, const char *typstr);

#endif
