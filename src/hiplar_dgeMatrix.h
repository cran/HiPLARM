#ifndef HIPLAR_GEMATRIX_H
#define HIPLAR_GEMATRIX_H

#include <R_ext/Boolean.h>
#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP hiplar_dgeMatrix_norm(SEXP obj, SEXP type);
SEXP hiplar_dgeMatrix_rcond(SEXP obj, SEXP type);

SEXP hiplar_dgeMatrix_crossprod(SEXP x, SEXP trans);
SEXP hiplar_dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP hiplar_dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP hiplar_dgeMatrix_LU(SEXP x, SEXP warn_singularity);
SEXP hiplar_dgeMatrix_determinant(SEXP x, SEXP logarithm);

SEXP hiplar_dgeMatrix_solve(SEXP a);
SEXP hiplar_dgeMatrix_matrix_solve(SEXP a, SEXP b);
SEXP hiplar_dgeMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);

#endif

