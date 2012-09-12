#include "plasma_dsyMatrix.h"
#include "P.h"
#include "plasma_init.h"
#include "hiplar_dbg.h"
#include <R_ext/Lapack.h>

double get_norm_sy(SEXP obj, const char *typstr)
{
#ifdef HIPLAR_WITH_PLASMA
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    double *work = (double *) NULL;

    typnm[0] = La_norm_type(typstr);

    if ((*typnm == 'F') && (CHECK_VVERSION_BEQ(2,4,5))) {
		error("not implemented");
    }

	if (*typnm == 'F') {
		work = (double *) R_alloc(2*R_PLASMA_NUM_THREADS, sizeof(double));
	} else {
		work = (double *) R_alloc(R_PLASMA_NUM_THREADS, sizeof(double));
    }

    return P_dlansy(typnm, uplo_P(obj),
        dims[0], REAL(GET_SLOT(obj, Matrix_xSym)), dims[0], work);

#endif
    return 0.0;
}

SEXP plasma_dsyMatrix_norm(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_PLASMA
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dsyMatrix_norm");
#endif

    return ScalarReal(get_norm_sy(obj, CHAR(asChar(type))));
#endif
    return R_NilValue;
}


SEXP plasma_dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rtP)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    int rt = asLogical(rtP); /* if(rt), compute b %*% a,  else  a %*% b */
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym)),
	m = bdims[0], n = bdims[1];
    double one = 1., zero = 0.;
    double *vx = REAL(GET_SLOT(val, Matrix_xSym));
    double *bcp = Memcpy(Alloca(m * n, double), vx, m * n);
    int info;
    R_CheckStack();

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dsyMatrix_matrix_mm");
#endif

    if ((rt && n != adims[0]) || (!rt && m != adims[0]))
	error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1) {
/* 	error(_("Matrices with zero extents cannot be multiplied")); */
    } else {

        info = P_dsymm(rt ? "R" :"L", uplo_P(a), m, n, one,
			REAL(GET_SLOT(a, Matrix_xSym)), adims[0], bcp,
			m, zero, vx, m);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dsymm", info);
        }
    }

    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}
