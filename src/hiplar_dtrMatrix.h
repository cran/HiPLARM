#ifndef HIPLAR_TRMATRIX_H
#define HIPLAR_TRMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP hiplar_dtrMatrix_solve(SEXP a);
SEXP hiplar_dtrMatrix_chol2inv(SEXP a);
SEXP hiplar_dtrMatrix_matrix_solve(SEXP a, SEXP b);
SEXP hiplar_dtrMatrix_matrix_mm   (SEXP a, SEXP b, SEXP right, SEXP trans);
SEXP hiplar_dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans);

#endif
