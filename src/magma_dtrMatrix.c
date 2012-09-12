#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#endif
#include "magma_dtrMatrix.h"
#include "setGPU.h"
#include "hiplar_dbg.h"

SEXP magma_dtrMatrix_solve(SEXP a)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = PROTECT(duplicate(a));
	int info, *Dim = INTEGER(GET_SLOT(val, Matrix_DimSym));
	double *A = REAL(GET_SLOT(val, Matrix_xSym)); 
	if(GPUFlag == 0) {
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtrMatrix Solve using dtrtri");
#endif
		F77_CALL(dtrtri)(uplo_P(val), diag_P(val), Dim,
			A , Dim, &info);
	}else if(GPUFlag == 1 && Interface == 0) {
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtrMatrix Solve using magma_dtrtri");
#endif
		magma_dtrtri(*uplo_P(val), *diag_P(val), Dim[0], A, Dim[0], &info);
	
	}else if(GPUFlag == 1 && Interface == 1) {
		double *d_A;
		cublasStatus retStatus;

#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtrMatrix Solve using magma_dtrtri_gpu");
#endif
		cublasAlloc(Dim[0] * Dim[1], sizeof(double), (void**)&d_A);
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector(Dim[0] * Dim[1], sizeof(double), A, 1, d_A, 1);

		magma_dtrtri_gpu(*uplo_P(val), *diag_P(val), Dim[0], d_A, Dim[0], &info);

		cublasGetVector(Dim[0] * Dim[1], sizeof(double), d_A, 1, A, 1);
		
		cublasFree(d_A);
	}else
		error(_("GPU Flag not set correctly"));

	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}


SEXP magma_dtrMatrix_chol2inv(SEXP a)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix")));
	int info, n;

	slot_dup(val, a, Matrix_DimSym);
	slot_dup(val, a, Matrix_uploSym);
	slot_dup(val, a, Matrix_diagSym);
	slot_dup(val, a, Matrix_DimNamesSym);
	slot_dup(val, a, Matrix_xSym);
	n = *INTEGER(GET_SLOT(val, Matrix_DimSym));
	double *A = REAL(GET_SLOT(val, Matrix_xSym));

	if(GPUFlag == 0) {
		F77_CALL(dpotri)(uplo_P(val), &n,
				A, &n, &info);
	}else if(GPUFlag == 1 && Interface == 0) {

#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtrMatrix chol2inv using magma_dpotri");
#endif
		magma_dpotri(*uplo_P(val), n, A, n, &info);

	}else if(GPUFlag == 1 && Interface == 1) {
		double *d_A;
		cublasStatus retStatus;

#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtrMatrix chol2inv using magma_dpotri_gpu");
#endif
		cublasAlloc(n * n, sizeof(double), (void**)&d_A);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation on Device"));
		/********************************************/

		cublasSetVector(n * n, sizeof(double), A, 1, d_A, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		magma_dpotri_gpu(*uplo_P(val), n, d_A, n, &info);

		cublasGetVector(n * n, sizeof(double), d_A, 1, A, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasFree(d_A);

	}else
		error(_("GPU Flag not set correctly"));


	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}


SEXP magma_dtrMatrix_matrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP ans = PROTECT(dup_mMatrix_as_dgeMatrix(b));
	int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	    *bdims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
	int n = bdims[0], nrhs = bdims[1];
	double one = 1.0;
	double *A =  REAL(GET_SLOT(a, Matrix_xSym));
	double *B = REAL(GET_SLOT(ans, Matrix_xSym));

	if (*adims != *bdims || bdims[1] < 1 || *adims < 1 || *adims != adims[1])
		error(_("Dimensions of system to be solved are inconsistent"));

	if(GPUFlag == 0) {
	
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtrMatrix, matrix solve using dtrsm");
#endif
		F77_CALL(dtrsm)("L", uplo_P(a), "N", diag_P(a),
				&n, &nrhs, &one, A, &n,
				B, &n);

	} else if(GPUFlag == 1) {
		double *d_A, *d_B;
		cublasStatus retStatus;

#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: dtrMatrix, matrix, solve using cublasDtrsm");
#endif
		cublasAlloc(n * n, sizeof(double), (void**)&d_A);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of A"));
		/********************************************/

		cublasAlloc(n * nrhs, sizeof(double), (void**)&d_B);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of B"));
		/********************************************/

		cublasSetVector(n * n, sizeof(double), A, 1, d_A, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector(n * nrhs, sizeof(double), B, 1, d_B, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		char U = *uplo_P(a);
		char D = *diag_P(a);

		cublasDtrsm('L', U, 'N', D, n, nrhs, one, d_A, n, d_B, nrhs);

		cublasGetVector(n * nrhs, sizeof(double), d_B, 1, B, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer From Device"));
		/********************************************/


		cublasFree(d_A);
		cublasFree(d_B);
	}else
		error(_("GPU Flag not set correctlt "));

	UNPROTECT(1);
	return ans;
#endif
	return R_NilValue;
}


SEXP magma_dtrMatrix_matrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans)
{
#ifdef HIPLAR_WITH_MAGMA
	/* Because a must be square, the size of the answer, val,
	 * is the same as the size of b */
	SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b));
	int rt = asLogical(right); /* if(rt), compute b %*% op(a),  else  op(a) %*% b */
	int tr = asLogical(trans);/* if true, use t(a) */
	int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	    *bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
	int m = bdims[0], n = bdims[1];
	double one = 1.;
	double *A = REAL(GET_SLOT(a, Matrix_xSym));
	double *B = REAL(GET_SLOT(val, Matrix_xSym));

	if (adims[0] != adims[1])
		error(_("dtrMatrix must be square"));
	if ((rt && adims[0] != n) || (!rt && adims[1] != m))
		error(_("Matrices are not conformable for multiplication"));
	if (m < 1 || n < 1) {
		error(_("Matrices with zero extents cannot be multiplied")); 
	} else if(GPUFlag == 0) {/* BLAS */
		F77_CALL(dtrmm)(rt ? "R" : "L", uplo_P(a),
				/*trans_A = */ tr ? "T" : "N",
				diag_P(a), &m, &n, &one,
				A , adims,
				B , &m);
	
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Performing dtrMatrix multiplication using dtrmm");
#endif
	} else if(GPUFlag == 1) {

		double *d_A, *d_B;
		cublasStatus retStatus;
		char uplo_a_ch_ = *uplo_P(a); //cublasDtrmm takes type char rather than const char *
		char diag_a_ch_ = *diag_P(a);

#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Performing dtrMatrix multiplication using cublasDtrmm");
#endif
		cublasAlloc(n * n, sizeof(double), (void**)&d_A);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of A"));
		/********************************************/

		cublasAlloc(m * n, sizeof(double), (void**)&d_B);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory allocation of B"));
		/********************************************/

		cublasSetVector(m * n, sizeof(double), B, 1, d_B, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector(n * n, sizeof(double), A, 1, d_A, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasDtrmm(rt ? 'R' : 'L', uplo_a_ch_, tr ? 'T' : 'N',  diag_a_ch_, m, n, one, d_A, adims[0], d_B, m);	

		cublasGetVector(m * n, sizeof(double), d_B, 1, B, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer from Device"));
		/********************************************/

		cublasFree(d_A);
		cublasFree(d_B);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error Freeing GPU memory"));
		/********************************************/
	
	}else
		error(_("Error with GPU Flag"));

	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}


SEXP magma_dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans)
{
#ifdef HIPLAR_WITH_MAGMA
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
	const char *uplo_a_ch = CHAR(STRING_ELT(uplo_a, 0)), /* = uplo_P(a) */
		*diag_a_ch = CHAR(STRING_ELT(diag_a, 0)); /* = diag_P(a) */

	Rboolean same_uplo = (*uplo_a_ch == *uplo_P(b)),
		   uDiag_b = /* -Wall: */ FALSE;
    int i;

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
		double alpha = 1.0;
		double *A = REAL(GET_SLOT(a, Matrix_xSym));
		double *B = REAL(GET_SLOT(val, Matrix_xSym));
		if(GPUFlag == 0) {
			/* Level 3 BLAS - DTRMM(): Compute one of the matrix multiplication operations
			 * B := alpha*op( A )*B ["L"], or B := alpha*B*op( A ) ["R"],
			 *	where trans_A determines  op(A):=  A   "N"one  or
			 *				  op(A):= t(A) "T"ransposed */
			F77_CALL(dtrmm)(rt ? "R" : "L", uplo_a_ch,
					/*trans_A = */ tr ? "T" : "N", diag_a_ch, &n, &n, &alpha,
					A , adims,
					B , &n);

#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Performing dtrMatrix multiplication using dtrmm");
#endif
		} else if(GPUFlag == 1) {

			double *d_A, *d_B;
			cublasStatus retStatus;
			char uplo_a_ch_ = *uplo_a_ch; //cublasDtrmm takes type char rather than const char *
			char diag_a_ch_ = *diag_a_ch;

#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Performing dtrMatrix multiplication using cublasDtrmm");
#endif
			cublasAlloc(n * n, sizeof(double), (void**)&d_A);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory allocation of A"));
			/********************************************/

			cublasAlloc(n * n, sizeof(double), (void**)&d_B);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory allocation of B"));
			/********************************************/

			cublasSetVector(n * n, sizeof(double), B, 1, d_B, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			cublasSetVector(n * n, sizeof(double), A, 1, d_A, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			cublasDtrmm(rt ? 'R' : 'L', uplo_a_ch_, tr ? 'T' : 'N',  diag_a_ch_, n, n, alpha, d_A, adims[0], d_B, n);	

			cublasGetVector(n * n, sizeof(double), d_B, 1, B, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer from Device"));
			/********************************************/

			cublasFree(d_A);
			cublasFree(d_B);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error Freeing GPU memory"));
			/********************************************/
		
		}
	}else	
		error(_("GPU Flag not set correctly"));

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
