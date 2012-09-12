#include "plasma_dtrMatrix.h"
#include "P.h"
#include "hiplar_dbg.h"

SEXP plasma_dtrMatrix_solve(SEXP a)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP val = PROTECT(duplicate(a));
    int info, *Dim = INTEGER(GET_SLOT(val, Matrix_DimSym));

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dtrMatrix_solve");
#endif

    info = P_dtrtri(uplo_P(val), diag_P(val),
        Dim[0], REAL(GET_SLOT(val, Matrix_xSym)), Dim[0]);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dtrtri", info);
    }

    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dtrMatrix_chol2inv(SEXP a)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix")));
    int info, n;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dtrMatrix_chol2inv");
#endif

    slot_dup(val, a, Matrix_DimSym);
    slot_dup(val, a, Matrix_uploSym);
    slot_dup(val, a, Matrix_diagSym);
    slot_dup(val, a, Matrix_DimNamesSym);
    slot_dup(val, a, Matrix_xSym);
    n = *INTEGER(GET_SLOT(val, Matrix_DimSym));

	info = P_dpotri(uplo_P(val), n, REAL(GET_SLOT(val, Matrix_xSym)), n);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dpotri", info);
    }

    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dtrMatrix_matrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP ans = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    int n = bdims[0], nrhs = bdims[1];
    double one = 1.0;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dtrMatrix_matrix_solve");
#endif

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1 || *adims != adims[1])
	error(_("Dimensions of system to be solved are inconsistent"));

	info = P_dtrsm("L", uplo_P(a), "N", diag_P(a), n, nrhs, one,
        REAL(GET_SLOT(a, Matrix_xSym)), n, REAL(GET_SLOT(ans, Matrix_xSym)), n);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dtrsm", info);
    }

    UNPROTECT(1);
    return ans;
#endif
    return R_NilValue;
}

/* to be used for all three: '%*%', crossprod() and tcrossprod() */
SEXP plasma_dtrMatrix_matrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans)
{
#ifdef HIPLAR_WITH_PLASMA
    /* Because a must be square, the size of the answer, val,
     * is the same as the size of b */
    SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    int rt = asLogical(right); /* if(rt), compute b %*% op(a),  else  op(a) %*% b */
    int tr = asLogical(trans);/* if true, use t(a) */
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int m = bdims[0], n = bdims[1];
    double one = 1.;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dtrMatrix_matrix_mm");
#endif

    if (adims[0] != adims[1])
	error(_("dtrMatrix must be square"));
    if ((rt && adims[0] != n) || (!rt && adims[1] != m))
	error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1) {
/* 	error(_("Matrices with zero extents cannot be multiplied")); */
	} else { /* BLAS */

        info = P_dtrmm(rt ? "R" : "L", uplo_P(a), tr ? "T" : "N",
			diag_P(a), m, n, one,
			REAL(GET_SLOT(a, Matrix_xSym)), adims[0],
			REAL(GET_SLOT(val, Matrix_xSym)), m);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dtrmm", info);
        }
    }

    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans)
{
#ifdef HIPLAR_WITH_PLASMA
    /* to be called from "%*%" and crossprod(), tcrossprod(),
     * from  ../R/products.R
     *
     * TWO cases : (1) result is triangular  <=> uplo are equal
     * ===         (2) result is "general"
     */
    SEXP val,/* = in case (2):  PROTECT(dup_mMatrix_as_dgeMatrix(b)); */
	d_a = GET_SLOT(a, Matrix_DimSym),
	uplo_a = GET_SLOT(a, Matrix_uploSym),
	diag_a = GET_SLOT(a, Matrix_diagSym);
    /* if(rt), compute b %*% a,  else  a %*% b */
    int rt = asLogical(right);
    int tr = asLogical(trans);/* if true, use t(a) */
    int *adims = INTEGER(d_a), n = adims[0];
    double *valx = (double *) NULL /*Wall*/;
    const char
	*uplo_a_ch = CHAR(STRING_ELT(uplo_a, 0)), /* = uplo_P(a) */
	*diag_a_ch = CHAR(STRING_ELT(diag_a, 0)); /* = diag_P(a) */
    Rboolean same_uplo = (*uplo_a_ch == *uplo_P(b)),
	uDiag_b = /* -Wall: */ FALSE;
    int i, info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dtrMatrix_dtrMatrix_mm");
#endif

    if (INTEGER(GET_SLOT(b, Matrix_DimSym))[0] != n)
	/* validity checking already "assures" square matrices ... */
	error(_("dtrMatrices in %*% must have matching (square) dim."));
    if(same_uplo) {
	/* ==> result is triangular -- "dtrMatrix" !
	 * val := dup_mMatrix_as_dtrMatrix(b) : */
	int sz = n * n;
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix")));
	SET_SLOT(val, Matrix_uploSym, duplicate(uplo_a));
	SET_SLOT(val, Matrix_DimSym,  duplicate(d_a));
	SET_DimNames(val, b);
	valx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz));
	Memcpy(valx, REAL(GET_SLOT(b, Matrix_xSym)), sz);
	if((uDiag_b = *diag_P(b) == 'U')) {
	    /* unit-diagonal b - may contain garbage in diagonal */
	    for (i = 0; i < n; i++)
		valx[i * (n+1)] = 1.;
	}
    } else { /* different "uplo" ==> result is "dgeMatrix" ! */
	val = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    }
    if (n >= 1) {
	double alpha = 1.;
	/* Level 3 BLAS - DTRMM(): Compute one of the matrix multiplication operations
	 * B := alpha*op( A )*B ["L"], or B := alpha*B*op( A ) ["R"],
	 *	where trans_A determines  op(A):=  A   "N"one  or
	 *				  op(A):= t(A) "T"ransposed */

	info = P_dtrmm(rt ? "R" : "L", uplo_a_ch,
			tr ? "T" : "N", diag_a_ch, n, n, alpha,
			REAL(GET_SLOT(a, Matrix_xSym)), adims[0],
			REAL(GET_SLOT(val, Matrix_xSym)), n);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dtrmm", info);
    }

    }
    if(same_uplo) {
	make_d_matrix_triangular(valx, a); /* set "other triangle" to 0 */
	if(*diag_a_ch == 'U' && uDiag_b) /* result remains uni-diagonal */
	    SET_SLOT(val, Matrix_diagSym, duplicate(diag_a));
    }
    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}
