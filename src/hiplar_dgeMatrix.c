#include "hiplar_dgeMatrix.h"
#include "plasma_dgeMatrix.h"
#include "magma_dgeMatrix.h"
#include "hiplar_at.h"
#include "hiplar_dbg.h"

SEXP hiplar_dgeMatrix_norm(SEXP obj, SEXP type) {
   	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_norm");
#endif

    size = dims[0]; /* dims[0] * dims[1] ? */

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dgeMatrix_norm))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_norm(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_norm(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_rcond(SEXP obj, SEXP type) {
   	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_rcond");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_rcond))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_rcond(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_rcond(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_crossprod(SEXP x, SEXP trans) {
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_crossprod");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dgeMatrix_crossprod))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_crossprod(x, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_crossprod(x, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans) {
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_dgeMatrix_crossprod");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_dgeMatrix_crossprod))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_dgeMatrix_crossprod(x, y, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_dgeMatrix_crossprod(x, y, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans) {
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_matrix_crossprod");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_matrix_crossprod))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_matrix_crossprod(x, y, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_matrix_crossprod(x, y, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_LU(SEXP x, SEXP warn_singularity) {
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_LU");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_LU))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_LU(x, warn_singularity);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_LU(x, warn_singularity);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_determinant(SEXP x, SEXP logarithm) {
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_determinant");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_determinant))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_determinant(x, logarithm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_determinant(x, logarithm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_solve(SEXP a) {
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_solve");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_solve))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_solve(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_solve(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_matrix_solve(SEXP a, SEXP b) {
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_matrix_solve");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_matrix_solve))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_matrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_matrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


SEXP hiplar_dgeMatrix_matrix_mm(SEXP a, SEXP b, SEXP right) {
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeMatrix_matrix_mm");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dgeMatrix_matrix_mm))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dgeMatrix_matrix_mm(a, b, right);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dgeMatrix_matrix_mm(a, b, right);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif
    return R_NilValue;

}


