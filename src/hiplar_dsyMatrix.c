#include "hiplar_dsyMatrix.h"
#include "plasma_dsyMatrix.h"
#include "magma_dsyMatrix.h"
#include "hiplar_at.h"
#include "hiplar_dbg.h"


SEXP hiplar_dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rt) {
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dsyMatrix_matrix_mm");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dsyMatrix_matrix_mm))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dsyMatrix_matrix_mm(a, b, rt);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dsyMatrix_matrix_mm(a, b, rt);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dsyMatrix_norm(SEXP obj, SEXP type) {
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    int size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dsyMatrix_norm");
#endif

    size = dims[0];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dsyMatrix_norm))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dsyMatrix_norm(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dsyMatrix_norm(obj, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif


    return R_NilValue;
}


