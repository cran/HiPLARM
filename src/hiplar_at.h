#ifndef HIPLAR_AT
#define HIPLAR_AT

#define HIPLAR_USE_PLASMA 1
#define HIPLAR_USE_MAGMA  2
#define HIPLAR_USE_AUTO   3


extern int hiplar_library;


/* dgeMatrix */

extern int xover_dgeMatrix_norm;
extern int xover_dgeMatrix_rcond;
extern int xover_dgeMatrix_crossprod;
extern int xover_dgeMatrix_dgeMatrix_crossprod;
extern int xover_dgeMatrix_matrix_crossprod;
extern int xover_dgeMatrix_LU;
extern int xover_dgeMatrix_determinant;
extern int xover_dgeMatrix_solve;
extern int xover_dgeMatrix_matrix_solve;
extern int xover_dgeMatrix_matrix_mm;


/* dpoMatrix */

extern int xover_dpoMatrix_rcond;
extern int xover_dpoMatrix_solve;
extern int xover_dpoMatrix_matrix_solve;
extern int xover_dpoMatrix_dgeMatrix_solve;
extern int xover_dpoMatrix_chol;


/* dsyMatrix */

extern int xover_dsyMatrix_matrix_mm;
extern int xover_dsyMatrix_norm;


/* dtrMatrix */

extern int xover_dtrMatrix_solve;
extern int xover_dtrMatrix_chol2inv;
extern int xover_dtrMatrix_matrix_solve;
extern int xover_dtrMatrix_matrix_mm;
extern int xover_dtrMatrix_dtrMatrix_mm;

#endif
