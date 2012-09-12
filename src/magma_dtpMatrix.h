#ifndef MAGMA_TPMATRIX_H
#define MAGMA_TPMATRIX_H

#include "Mutils.h"
#include <R_ext/Lapack.h>

SEXP magma_dtpMatrix_matrix_mm(SEXP x, SEXP y);
SEXP magma_dgeMatrix_dtpMatrix_mm(SEXP x, SEXP y);
SEXP magma_dtpMatrix_matrix_solve(SEXP a, SEXP b);

#endif
