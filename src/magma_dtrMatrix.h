#ifndef MAGMA_TRMATRIX_H
#define MAGMA_TRMATRIX_H

#include "Mutils.h"
#include <R_ext/Lapack.h>

SEXP magma_dtrMatrix_solve(SEXP a);
SEXP magma_dtrMatrix_chol2inv(SEXP a);
SEXP magma_dtrMatrix_matrix_solve(SEXP a, SEXP b);
SEXP magma_dtrMatrix_matrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans);
SEXP magma_dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans);

#endif
