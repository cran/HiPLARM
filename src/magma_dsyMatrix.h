#ifndef MAGMA_SYMATRIX_H
#define MAGMA_SYMATRIX_H

#include "Mutils.h"
#include <R_ext/Lapack.h>

SEXP magma_dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rtP);
SEXP magma_dsyMatrix_norm(SEXP obj, SEXP type);

#endif
