#ifndef PLASMA_GEMATRIX_H
#define PLASMA_GEMATRIX_H

#include <R_ext/Boolean.h>
#include <R_ext/Lapack.h>
#include "Mutils.h"


SEXP plasma_dgeMatrix_norm(SEXP obj, SEXP norm);
SEXP plasma_dgeMatrix_rcond(SEXP obj, SEXP type);

/* for crossprod() and tcrossprod() : */
SEXP plasma_dgeMatrix_crossprod(SEXP x, SEXP trans);
SEXP plasma_dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP plasma_dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans);

SEXP plasma_dgeMatrix_LU (SEXP x, SEXP warn_singularity);
SEXP plasma_dgeMatrix_LU_(SEXP x, Rboolean warn_sing);
SEXP plasma_dgeMatrix_determinant(SEXP x, SEXP logarithm);

SEXP dgeMatrix_Schur(SEXP x, SEXP vectors);

SEXP plasma_dgeMatrix_solve(SEXP a);
SEXP plasma_dgeMatrix_matrix_solve(SEXP a, SEXP b);
SEXP plasma_dgeMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);

SEXP dgeMatrix_svd(SEXP x, SEXP nu, SEXP nv);
SEXP dgeMatrix_exp(SEXP x);
SEXP dgeMatrix_colsums(SEXP x, SEXP naRmP, SEXP cols, SEXP mean);

/* DGESDD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or      */
/* right singular vectors.  If singular vectors are desired, it uses a */
/* divide-and-conquer algorithm.                                   */
void F77_NAME(dgesdd)(const char *jobz,
		      const int *m, const int *n,
		      double *a, const int *lda, double *s,
		      double *u, const int *ldu,
		      double *vt, const int *ldvt,
		      double *work, const int *lwork, int *iwork, int *info);


#endif

