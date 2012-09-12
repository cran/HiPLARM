#ifndef HIPLAR_SYMATRIX_H
#define HIPLAR_SYMATRIX_H

#include "Mutils.h"
#include "R_ext/Lapack.h"

SEXP hiplar_dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rt);
SEXP hiplar_dsyMatrix_norm(SEXP obj, SEXP type);

#endif
