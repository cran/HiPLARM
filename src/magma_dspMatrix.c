#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#endif

#include "magma_dspMatrix.h"
#include "setGPU.h"
#include "hiplar_dbg.h"

SEXP magma_dspMatrix_matrix_mm(SEXP a, SEXP b)
{ 
#ifdef HIPLAR_WITH_MAGMA

	SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b));
	int *bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
	int i, ione = 1, n = bdims[0], nrhs = bdims[1];
	const char *uplo = uplo_P(a);
	double *ax = REAL(GET_SLOT(a, Matrix_xSym)), one = 1., zero = 0.,
		 *vx = REAL(GET_SLOT(val, Matrix_xSym));
	double *bx = Alloca(n * nrhs, double);
	R_CheckStack();

	Memcpy(bx, vx, n * nrhs);
	if (bdims[0] != n)
		error(_("Matrices are not conformable for multiplication"));
	if (nrhs < 1 || n < 1) {
		//	This was commented out
			error(_("Matrices with zero extents cannot be multiplied")); 
	} else if(GPUFlag == 0) {	
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Performing cross prod with dspmv");
#endif
		for (i = 0; i < nrhs; i++)	
			F77_CALL(dspmv)(uplo, &n, &one, ax, bx + i * n, &ione,
					&zero, vx + i * n, &ione);
	}
	else if(GPUFlag == 1) {
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Performing cross prod with cublasDspmv");
#endif
		double *d_ax, *d_bx, *d_vx;
		cublasStatus retStatus;
		cublasAlloc(n * (n + 1) / 2, sizeof(double), (void**)&d_ax);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation"));
		/********************************************/

		cublasAlloc(n * nrhs, sizeof(double), (void**)&d_bx);	

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation"));
		/********************************************/

		cublasAlloc(n * nrhs, sizeof(double), (void**)&d_vx);	

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation"));
		/********************************************/

		cublasSetVector( n * (n + 1) / 2 , sizeof(double), ax, 1, d_ax, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector( n * nrhs, sizeof(double), bx, 1, d_bx, 1 );

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector( n * nrhs, sizeof(double), vx, 1, d_vx, 1 );

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		//for (i = 0; i < nrhs; i++)	
		//F77_CALL(dspmv)(uplo, &n, &one, ax, bx + i * n, &ione,
		//	&zero, vx + i * n, &ione);

		char U = *uplo_P(a);
		for(i = 0; i < nrhs; i++)
			cublasDspmv(U, n, one, d_ax, d_bx + i * n, ione, zero, d_vx + i * n, ione);

		cublasGetVector( n * nrhs, sizeof(double), d_vx, 1, vx, 1 );

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasFree(d_ax);
		cublasFree(d_bx);
		cublasFree(d_vx);

	}else
		error(_("Error with GPU Flag"));

	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}
