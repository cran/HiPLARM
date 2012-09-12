#ifndef PLASMA_SYMATRIX_H
#define PLASMA_SYMATRIX_H

#include "Mutils.h"
#include "R_ext/Lapack.h"

SEXP plasma_dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rt);
SEXP plasma_dsyMatrix_norm(SEXP obj, SEXP type);
double get_norm_sy(SEXP obj, const char *typstr);

#endif
