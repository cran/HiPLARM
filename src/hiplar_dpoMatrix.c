#include "hiplar_dpoMatrix.h"
#include "plasma_dpoMatrix.h"
#include "magma_dpoMatrix.h"
#include "hiplar_at.h"
#include "hiplar_dbg.h"


SEXP hiplar_dpoMatrix_rcond(SEXP obj, SEXP type) {
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dpoMatrix_rcond");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if ((hiplar_library == HIPLAR_USE_PLASMA) ||
        ((hiplar_library == HIPLAR_USE_AUTO) &&  (size < xover_dpoMatrix_rcond))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dpoMatrix_rcond(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dpoMatrix_rcond(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dpoMatrix_solve(SEXP a) {
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dpoMatrix_solve");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dpoMatrix_solve))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dpoMatrix_solve(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dpoMatrix_solve(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dpoMatrix_matrix_solve(SEXP a, SEXP b) {
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dpoMatrix_matrix_solve");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dpoMatrix_matrix_solve))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dpoMatrix_matrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dpoMatrix_matrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dpoMatrix_dgeMatrix_solve(SEXP a, SEXP b) {
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dpoMatrix_dgeMatrix_solve");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dpoMatrix_dgeMatrix_solve))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dpoMatrix_dgeMatrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dpoMatrix_dgeMatrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dpoMatrix_chol(SEXP x) {
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dpoMatrix_chol");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dpoMatrix_chol))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dpoMatrix_chol(x);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dpoMatrix_chol(x);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


