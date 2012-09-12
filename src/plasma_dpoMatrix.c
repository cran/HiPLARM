#include "plasma_dpoMatrix.h"
#include "P.h"
#include "hiplar_dbg.h"


SEXP plasma_dpoMatrix_chol(SEXP x)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP val = get_factors(x, "Cholesky"),
	dimP = GET_SLOT(x, Matrix_DimSym),
	uploP = GET_SLOT(x, Matrix_uploSym);
    const char *uplo = CHAR(STRING_ELT(uploP, 0));
    int *dims = INTEGER(dimP), info;
    int n = dims[0];
    double *vx;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dpoMatrix_chol");
#endif

    if (val != R_NilValue) {
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: plasma_dpoMatrix_chol using cached value");
#endif
        return val;
    }

    dims = INTEGER(dimP);
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("Cholesky")));
    SET_SLOT(val, Matrix_uploSym, duplicate(uploP));
    SET_SLOT(val, Matrix_diagSym, mkString("N"));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n * n));
    AZERO(vx, n * n);
    F77_CALL(dlacpy)(uplo, &n, &n, REAL(GET_SLOT(x, Matrix_xSym)), &n, vx, &n);
    if (n > 0) {
		info = P_dpotrf(uplo, n, vx, n);
    	if (info) {
    	    if(info > 0) {
        		error(_("the leading minor of order %d is not positive definite"), info);
            } else {
        		error(_("PLASMA routine %s returned error code %d"), "PLASMA_dpotrf", info);
            }
	    }

    }
    UNPROTECT(1);
    return set_factors(x, val, "Cholesky");
#endif
    return R_NilValue;
}

SEXP plasma_dpoMatrix_rcond(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP Chol = plasma_dpoMatrix_chol(obj);
    const char typnm[] = {'O', '\0'};	/* always use the one norm */
    int *dims = INTEGER(GET_SLOT(Chol, Matrix_DimSym)), info;
    double anorm = get_norm_sy(obj, typnm), rcond;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dpoMatrix_rcond");
#endif

    F77_CALL(dpocon)(uplo_P(Chol),
		     dims, REAL(GET_SLOT(Chol, Matrix_xSym)),
		     dims, &anorm, &rcond,
		     (double *) R_alloc(3*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info);
    return ScalarReal(rcond);
#endif
    return R_NilValue;
}

SEXP plasma_dpoMatrix_solve(SEXP x)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP Chol = plasma_dpoMatrix_chol(x);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix")));
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dpoMatrix_solve");
#endif

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    slot_dup(val, Chol, Matrix_uploSym);
    slot_dup(val, Chol, Matrix_xSym);
    slot_dup(val, Chol, Matrix_DimSym);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(x, Matrix_DimNamesSym)));

	info = P_dpotri(uplo_P(val), dims[0], REAL(GET_SLOT(val, Matrix_xSym)), dims[0]);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dpotri", info);
    }

    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dpoMatrix_dgeMatrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP Chol = plasma_dpoMatrix_chol(a),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dpoMatrix_dgeMatrix_solve");
#endif

    if (adims[1] != bdims[0])
	error(_("Dimensions of system to be solved are inconsistent"));
    if (adims[0] < 1 || bdims[1] < 1)
	error(_("Cannot solve() for matrices with zero extents"));
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    slot_dup(val, b, Matrix_DimSym);
    slot_dup(val, b, Matrix_xSym);

    info = P_dpotrs(uplo_P(Chol), adims[0], bdims[1],
        REAL(GET_SLOT(Chol, Matrix_xSym)), adims[0],
        REAL(GET_SLOT(val, Matrix_xSym)), bdims[0]);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dpotrs", info);
    }

    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dpoMatrix_matrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP Chol = plasma_dpoMatrix_chol(a),
	val = PROTECT(duplicate(b));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(getAttrib(b, R_DimSymbol)),
	info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dpoMatrix_matrix_solve");
#endif

    if (!(isReal(b) && isMatrix(b)))
	error(_("Argument b must be a numeric matrix"));
    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));

    info = P_dpotrs(uplo_P(Chol), adims[0], bdims[1],
        REAL(GET_SLOT(Chol, Matrix_xSym)), adims[0],
        REAL(val), bdims[0]);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dpotrs", info);
    }

    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}
