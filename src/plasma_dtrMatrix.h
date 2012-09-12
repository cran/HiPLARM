#ifndef PLASMA_TRMATRIX_H
#define PLASMA_TRMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP plasma_dtrMatrix_solve(SEXP a);
SEXP plasma_dtrMatrix_chol2inv(SEXP a);
SEXP plasma_dtrMatrix_matrix_solve(SEXP a, SEXP b);
SEXP plasma_dtrMatrix_matrix_mm   (SEXP a, SEXP b, SEXP right, SEXP trans);
SEXP plasma_dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans);

#endif
