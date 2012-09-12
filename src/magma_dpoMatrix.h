#ifndef MAGMA_POMATRIX_H
#define MAGMA_POMATRIX_H

#include "Mutils.h"
#include <R_ext/Lapack.h>

SEXP magma_dpoMatrix_rcond(SEXP obj, SEXP type);
SEXP magma_dpoMatrix_solve(SEXP x);
SEXP magma_dpoMatrix_matrix_solve(SEXP a, SEXP b);
SEXP magma_dpoMatrix_dgeMatrix_solve(SEXP a, SEXP b);
SEXP magma_dpoMatrix_chol(SEXP x);
double magma_get_norm_sy(SEXP obj, const char *typstr);
#endif

