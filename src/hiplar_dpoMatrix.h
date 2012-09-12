#ifndef HIPLAR_POMATRIX_H
#define HIPLAR_POMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP hiplar_dpoMatrix_rcond(SEXP obj, SEXP type);
SEXP hiplar_dpoMatrix_solve(SEXP a);
SEXP hiplar_dpoMatrix_matrix_solve(SEXP a, SEXP b);
SEXP hiplar_dpoMatrix_dgeMatrix_solve(SEXP a, SEXP b);
SEXP hiplar_dpoMatrix_chol(SEXP x);

#endif
