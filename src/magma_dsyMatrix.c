#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#endif
#include "magma_dsyMatrix.h"
#include "magma_dpoMatrix.h"
#include "setGPU.h"
#include "hiplar_dbg.h"


SEXP magma_dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rtP)
{
#ifdef HIPLAR_WITH_MAGMA

    SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    int rt = asLogical(rtP); /* if(rt), compute b %*% a,  else  a %*% b */
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
        *bdims = INTEGER(GET_SLOT(val, Matrix_DimSym)),
        m = bdims[0], n = bdims[1];
    double one = 1., zero = 0.;
    double *vx = REAL(GET_SLOT(val, Matrix_xSym));
    double *bcp = Memcpy(Alloca(m * n, double), vx, m * n);

    double *A = REAL(GET_SLOT(a, Matrix_xSym));

    R_CheckStack();
    
		if ((rt && n != adims[0]) || (!rt && m != adims[0]))
        error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1) {
        //	This was commented out
        //    error(_("Matrices with zero extents cannot be multiplied")); 
    } else if(GPUFlag == 1) {
        double *d_A, *d_B, *d_C;
        cublasStatus retStatus;
        int i;
#ifdef HIPLAR_DBG
        R_ShowMessage("DBG: Performing dsy||dpo crossprod using magmablas_dsymm");
#endif
        cublasAlloc(adims[0] * adims[1], sizeof(double), (void**)&d_A);

        /* Error Checking */
        retStatus = cublasGetError ();
        if (retStatus != CUBLAS_STATUS_SUCCESS) 
            error(_("CUBLAS: Error in Memory Allocation"));
        /********************************************/

        cublasAlloc(m * n, sizeof(double), (void**)&d_B);	

        /* Error Checking */
        retStatus = cublasGetError ();
        if (retStatus != CUBLAS_STATUS_SUCCESS) 
            error(_("CUBLAS: Error in Memory Allocation"));
        /********************************************/

        cublasAlloc(m * n, sizeof(double), (void**)&d_C);	

        /* Error Checking */
        retStatus = cublasGetError ();
        if (retStatus != CUBLAS_STATUS_SUCCESS) 
            error(_("CUBLAS: Error in Memory Allocation"));
        /********************************************/

        cublasSetVector( adims[0]  * adims[1] , sizeof(double), A, 1, d_A, 1);

        /* Error Checking */
        retStatus = cublasGetError ();
        if (retStatus != CUBLAS_STATUS_SUCCESS) 
            error(_("CUBLAS: Error in Data Transfer to Device"));
        /********************************************/

        cublasSetVector( m * n, sizeof(double), bcp, 1, d_B, 1 );

        /* Error Checking */
        retStatus = cublasGetError ();
        if (retStatus != CUBLAS_STATUS_SUCCESS) 
            error(_("CUBLAS: Error in Data Transfer to Device"));
        /********************************************/

        char U = *uplo_P(a);
        cublasDsymm(rt ? 'R' : 'L', U, m, n, one, d_A, adims[0], d_B, m, zero, d_C,  m);
        cublasGetVector( m * n , sizeof(double), d_C, 1, vx, 1);

        /* Error Checking */
        retStatus = cublasGetError ();
        if (retStatus != CUBLAS_STATUS_SUCCESS) 
            error(_("CUBLAS: Error in Data Transfer from Device"));
        /********************************************/

        cublasFree(d_A);
        cublasFree(d_B);
        cublasFree(d_C);
    } else {
        F77_CALL(dsymm)(rt ? "R" :"L", uplo_P(a), &m, &n, &one,
                A , adims, bcp,
                &m, &zero, vx, &m);
        
#ifdef HIPLAR_DBG
        R_ShowMessage("DBG: Performing dsy||dpo crossprod using dsymm");
#endif
    }
    UNPROTECT(1);
    return val;
#endif
	return R_NilValue;
}

SEXP magma_dsyMatrix_norm(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_MAGMA
	      return ScalarReal(magma_get_norm_sy(obj, CHAR(asChar(type))));
#endif
	return R_NilValue;
}

