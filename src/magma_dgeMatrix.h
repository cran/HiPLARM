#ifndef MAGMA_GEMATRIX_H
#define MAGMA_GEMATRIX_H

#include "Mutils.h"
#include <R_ext/Lapack.h>

SEXP magma_dgeMatrix_norm(SEXP obj, SEXP type);
SEXP magma_dgeMatrix_rcond(SEXP obj, SEXP type);

SEXP magma_dgeMatrix_crossprod(SEXP x, SEXP trans);
SEXP magma_dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP magma_dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP magma_dgeMatrix_LU_(SEXP x, Rboolean warn_sing);
SEXP magma_dgeMatrix_LU(SEXP x, SEXP warn_singularity);

SEXP magma_dgeMatrix_solve(SEXP a);
SEXP magma_dgeMatrix_matrix_solve(SEXP a, SEXP b);
SEXP magma_dgeMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);

SEXP magma_dgeMatrix_determinant(SEXP x, SEXP logarithm);

#endif
