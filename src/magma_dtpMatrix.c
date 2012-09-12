#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#endif
#include "magma_dtpMatrix.h"
#include "setGPU.h"
#include "hiplar_dbg.h"

SEXP magma_dtpMatrix_matrix_mm(SEXP x, SEXP y)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(y));
	/* Since 'x' is square (n x n ),   dim(x %*% y) = dim(y) */
	int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	    *yDim = INTEGER(GET_SLOT(val, Matrix_DimSym));
	int ione = 1, j;
	const char *uplo = uplo_P(x), *diag = diag_P(x);
	double *xx = REAL(GET_SLOT(x, Matrix_xSym)),
		 *vx = REAL(GET_SLOT(val, Matrix_xSym));

	if (yDim[0] != xDim[1])
		error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
				xDim[0], xDim[1], yDim[0], yDim[1]);

	if(GPUFlag == 0) {	
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtpMatrix_matrix * matrix using dtpmv");
#endif
		for (j = 0; j < yDim[1]; j++) /* X %*% y[,j]  via BLAS 2 DTPMV(.) */
			F77_CALL(dtpmv)(uplo, "N", diag, yDim, xx,
					vx + j * yDim[0], &ione);
	}
	else if(GPUFlag == 1) {

		double *d_xx, *d_vx;	
		cublasStatus retStatus;
#ifdef HIPLAR_DBG		
		R_ShowMessage("DBG: dtpMatrix_matrix * matrix using cublasDtpmv");
#endif		
		cublasAlloc(xDim[0] * (xDim[1] + 1) / 2, sizeof(double), (void**)&d_xx);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of A"));
		/********************************************/

		cublasAlloc(yDim[0] * yDim[1], sizeof(double), (void**)&d_vx);
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of B"));
		/********************************************/

		//Ok to here	
		int n = xDim[0];
		cublasSetVector( n * (n + 1) / 2, sizeof(double), xx, 1, d_xx, 1);
		//cublasSetVector( xDim[1], sizeof(double), xx, 1, d_xx, 1);
		
		//cublasSetMatrix(xDim[0] - 1, xDim[1] - 1, sizeof(double), xx, xDim[0], d_xx, xDim[0]);	
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in %s Data Transfer to Device"), cudaGetErrorString(retStatus));
		/********************************************/

		cublasSetVector(yDim[0] * yDim[1], sizeof(double), vx, 1, d_vx, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in %s Data Transfer to Device"), cudaGetErrorString(retStatus));
		/********************************************/

		char U = *uplo_P(x);
		char D = *diag_P(x);

		for (j = 0; j < yDim[1]; j++)
			cublasDtpmv(U, 'N', D, yDim[0], d_xx, d_vx + j * yDim[0], ione);

		cublasGetVector(yDim[0] * yDim[1], sizeof(double), d_vx, 1, vx, 1);
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer from Device"));
		/********************************************/


		cublasFree(d_xx);
		cublasFree(d_vx);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Freeing"));
		/********************************************/
		//cublasDestroy(handle);

	}else
		error(_("GPU Flag not set correctly"));

	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}


SEXP magma_dgeMatrix_dtpMatrix_mm(SEXP x, SEXP y)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = PROTECT(duplicate(x));
	/* Since 'y' is square (n x n ),   dim(x %*% y) = dim(x) */
	int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	    *yDim = INTEGER(GET_SLOT(y, Matrix_DimSym));
	int i;
	const char *uplo = uplo_P(y), *diag = diag_P(y);
	double *yx = REAL(GET_SLOT(y, Matrix_xSym)),
		 *vx = REAL(GET_SLOT(val, Matrix_xSym));

	if (yDim[0] != xDim[1])
		error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
				xDim[0], xDim[1], yDim[0], yDim[1]);
	if(GPUFlag == 0) {
		
		for (i = 0; i < xDim[0]; i++)/* val[i,] := Y' %*% x[i,]  */
			F77_CALL(dtpmv)(uplo, "T", diag, yDim, yx,
					vx + i, /* incr = */ xDim);

#ifdef HIPLAR_DBG		
		R_ShowMessage("DBG: Performing dtrMatrix multiplication using dtpmv");
#endif
	}else if(GPUFlag == 1) {

		double *d_yx, *d_vx;	
		cublasStatus retStatus;
	
#ifdef HIPLAR_DBG		
		R_ShowMessage("DBG: dtpMatrix_matrix * dgeMatrix using cublasDtpmv");
#endif
		cublasAlloc(yDim[0] * yDim[1], sizeof(double), (void**)&d_yx);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of A"));
		/********************************************/

		cublasAlloc(xDim[0] * xDim[1], sizeof(double), (void**)&d_vx);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of B"));
		/********************************************/

		cublasSetVector(yDim[0] * yDim[1], sizeof(double), yx, 1, d_yx, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector(xDim[0] * xDim[1], sizeof(double), vx, 1, d_vx, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		char U = *uplo_P(y);
		char D = *diag_P(y);
		for( i = 0; i < xDim[0]; i++)
		cublasDtpmv(U, 'T', D, yDim[0], d_yx, d_vx + i, xDim[0]);
		
		cublasGetVector(xDim[0] * xDim[1], sizeof(double), d_vx, 1, vx, 1);

		cublasFree(d_yx);
		cublasFree(d_vx);

	}

	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}


SEXP magma_dtpMatrix_matrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b));
	/* Since 'x' is square (n x n ),   dim(x %*% b) = dim(b) */
	int *aDim = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	    *bDim = INTEGER(GET_SLOT(val, Matrix_DimSym));
	int ione = 1, j;
	const char *uplo = uplo_P(a), *diag = diag_P(a);
	double *ax = REAL(GET_SLOT(a, Matrix_xSym)),
		 *vx = REAL(GET_SLOT(val, Matrix_xSym));

	if (bDim[0] != aDim[1])
		error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
				aDim[0], aDim[1], bDim[0], bDim[1]);
	if(GPUFlag == 0) {
		for (j = 0; j < bDim[1]; j++) /* a^{-1} %*% b[,j]  via BLAS 2 DTPSV(.) */
			F77_CALL(dtpsv)(uplo, "N", diag, bDim, ax,
					vx + j * bDim[0], &ione);
	}else if(GPUFlag == 1) {
		double *d_ax, *d_vx;	
		cublasStatus retStatus;
	
#ifdef HIPLAR_DBG		
		R_ShowMessage("DBG: dtpMatrix_matrix * matrix using cublasDtpmv");
#endif		
		cublasAlloc(aDim[0] * aDim[1], sizeof(double), (void**)&d_ax);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of A"));
		/********************************************/

		cublasAlloc(bDim[0] * bDim[1], sizeof(double), (void**)&d_vx);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of B"));
		/********************************************/

		cublasSetVector(aDim[0] * aDim[1], sizeof(double), ax, 1, d_ax, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector(bDim[0] * bDim[1], sizeof(double), vx, 1, d_vx, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		char U = *uplo_P(a);
		char D = *diag_P(a);

		cublasDtpsv(U, 'N', D, bDim[0], d_ax, d_vx, ione);

		cublasGetVector(bDim[0] * bDim[1], sizeof(double), d_vx, 1, vx, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer From Device"));
		/********************************************/


		cublasFree(d_ax);
		cublasFree(d_vx);
	}

	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}
