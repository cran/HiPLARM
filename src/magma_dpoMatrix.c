#include "magma_dpoMatrix.h"
#include "setGPU.h"
#include "hiplar_dbg.h"

#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#include <cublas.h>
#endif


/* Notes *
 *
 *
 * uplo_P(_x_) defined as CHAR(STRING_ELT(GET_SLOT(_x_, Matrix_uploSym), 0))
 */


double magma_get_norm_sy(SEXP obj, const char *typstr)
{
#ifdef HIPLAR_WITH_MAGMA
	char typnm[] = {'\0', '\0'};
	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
	double *work = (double *) NULL;
	int N = dims[0];
	int lda = N;
	double *A = REAL(GET_SLOT(obj, Matrix_xSym));
	typnm[0] = La_norm_type(typstr);

	const char *c = uplo_P(obj);

	//Magmablas dlansy only does I & M norms
	if(GPUFlag == 1 && (*typnm == 'I' || *typnm == 'M')) {
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Performing norm using magmablas_dlansy"); 
#endif
		double *dwork, *d_A, maxnorm;
		cublasAlloc(N, sizeof(double), (void**)&dwork);
		cublasAlloc(lda * N, sizeof(double), (void**)&d_A);
		cublasSetVector(N * lda, sizeof(double), A, 1, d_A, 1);
		maxnorm = magmablas_dlansy(typnm[0], *c ,N, d_A, lda, dwork);
		cublasFree(d_A);
		cublasFree(dwork);
		return maxnorm;
	}
	else {

		if (*typnm == 'I' || *typnm == 'O') {
			work = (double *) R_alloc(dims[0], sizeof(double));
		}

		return F77_CALL(dlansy)(typnm, uplo_P(obj),
				dims, A,
				dims, work);
	}
#endif
	return 0.0;
}


SEXP magma_dpoMatrix_chol(SEXP x)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = get_factors(x, "Cholesky"),
			 dimP = GET_SLOT(x, Matrix_DimSym),
			 uploP = GET_SLOT(x, Matrix_uploSym);

	const char *uplo = CHAR(STRING_ELT(uploP, 0));
	int *dims = INTEGER(dimP), info;
	int n = dims[0];
	double *vx;
	if (val != R_NilValue) return val;
	dims = INTEGER(dimP);
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("Cholesky")));
	SET_SLOT(val, Matrix_uploSym, duplicate(uploP));
	SET_SLOT(val, Matrix_diagSym, mkString("N"));
	SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
	vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n * n));
	AZERO(vx, n * n);
	
	//we could put in magmablas_dlacpy but it only
	//copies all of the matrix 
	F77_CALL(dlacpy)(uplo, &n, &n, REAL(GET_SLOT(x, Matrix_xSym)), &n, vx, &n);
	if (n > 0) {

		if(GPUFlag == 0){
#ifdef HIPLAR_DBG	
		R_ShowMessage("DBG: Cholesky decomposition using dpotrf;");
#endif
			F77_CALL(dpotrf)(uplo, &n, vx, &n, &info);
		}
		else if(GPUFlag == 1 && Interface == 0){
		
#ifdef HIPLAR_DBG	
			R_ShowMessage("DBG: Cholesky decomposition using magma_dpotrf;");
#endif			
			double *d_c;
			int nrows, ncols;
			nrows = ncols = n;
			int N2 = nrows * ncols;

			magma_int_t lda;
			lda = nrows;

			magma_dpotrf(uplo[0], ncols, vx, lda, &info);
		}
		else if(GPUFlag == 1 && Interface == 1) {
	
#ifdef HIPLAR_DBG	
			R_ShowMessage("DBG: Cholesky decomposition using magma_dpotrf_gpu;");	
#endif
			double *d_c;
			int nrows, ncols;
			nrows = ncols = n;
			int N2 = nrows * ncols;
			cublasStatus retStatus;

			magma_int_t lda;
			lda = nrows;

			cublasAlloc(lda * ncols, sizeof(double), (void**)&d_c);
			
			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory Allocation"));
			/********************************************/

			cublasSetVector(N2, sizeof(double), vx, 1, d_c, 1);
			
			/* Error Checking */
			retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Date Transfer to Device"));
			/********************************************/


			magma_dpotrf_gpu(uplo[0], ncols, d_c, lda, &info);
			
			cublasGetVector(nrows * ncols, sizeof(double), d_c, 1, vx, 1);		
			
			/* Error Checking */
			retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Date Transfer from Device"));
			/********************************************/
			
			cublasFree(d_c);
		}
		else
			error(_("MAGMA/LAPACK/Interface Flag not defined correctly"));
		}
		
	if (info) {
			if(info > 0)
				error(_("the leading minor of order %d is not positive definite"),
						info);
			else /* should never happen! */
				error(_("Lapack routine %s returned error code %d"), "dpotrf", info);
		}

	UNPROTECT(1);
	return set_factors(x, val, "Cholesky");
#endif
	return R_NilValue;
}

//	need to implement magma here, but as noted below
//	we are limited to the type of norm we can use (see get_norm_sy)
SEXP magma_dpoMatrix_rcond(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_MAGMA
    SEXP Chol = magma_dpoMatrix_chol(obj);
    const char typnm[] = {'O', '\0'};	// always use the one norm 
    int *dims = INTEGER(GET_SLOT(Chol, Matrix_DimSym)), info;
    double anorm = magma_get_norm_sy(obj, typnm), rcond;

    F77_CALL(dpocon)(uplo_P(Chol),
		     dims, REAL(GET_SLOT(Chol, Matrix_xSym)),
		     dims, &anorm, &rcond,
		     (double *) R_alloc(3*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info);
    return ScalarReal(rcond);
#endif
	return R_NilValue;
}

SEXP magma_dpoMatrix_solve(SEXP x)
{
#ifdef HIPLAR_WITH_MAGMA
    SEXP Chol = magma_dpoMatrix_chol(x);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix")));
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), info;

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    slot_dup(val, Chol, Matrix_uploSym);
    slot_dup(val, Chol, Matrix_xSym);
    slot_dup(val, Chol, Matrix_DimSym);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(x, Matrix_DimNamesSym)));
    double *A = REAL(GET_SLOT(val, Matrix_xSym));
		int N = *dims;	
		int lda = N;
		const char *uplo = uplo_P(val);
		
		if(GPUFlag == 0) {
			
			F77_CALL(dpotri)(uplo_P(val), dims, A, dims, &info);
		
		}
		else if(GPUFlag == 1 && Interface == 0) {
#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Solving using magma_dpotri");
#endif
			magma_dpotri(uplo[0], N, A, lda, &info);
		}
		else if(GPUFlag == 1 && Interface == 1){
			double *d_A;
			cublasStatus retStatus;
			cublasAlloc( N * lda , sizeof(double), (void**)&d_A);
#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Solving using magma_dpotri_gpu");
#endif		
			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory Allocation"));
			/********************************************/

			cublasSetVector( N * lda, sizeof(double), A, 1, d_A, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			magma_dpotri_gpu(uplo[0], N, d_A, lda, &info);

			cublasGetVector(N * lda, sizeof(double), d_A, 1, val, 1);
			
			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer from Device"));
			/********************************************/
			
			cublasFree(d_A);
		}
		else
			error(_("MAGMA/LAPACK/Interface Flag not defined correctly"));
	
		if (info) {
			if(info > 0)
				error(_("the leading minor of order %d is not positive definite"),
						info);
			else /* should never happen! */
				error(_("Lapack routine %s returned error code %d"), "dpotrf", info);
		}
		
		UNPROTECT(1);
    return val;
#endif
	return R_NilValue;
}


SEXP magma_dpoMatrix_dgeMatrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP Chol = magma_dpoMatrix_chol(a),
			 val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
	int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
			*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
			info;

	/* Checking Matrix Dimensions */
	if (adims[1] != bdims[0])
		error(_("Dimensions of system to be solved are inconsistent"));
	if (adims[0] < 1 || bdims[1] < 1)
		error(_("Cannot solve() for matrices with zero extents"));
	/* ****************************************** */
	
	SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
	slot_dup(val, b, Matrix_DimSym);
	slot_dup(val, b, Matrix_xSym);

	double *A = REAL(GET_SLOT(Chol, Matrix_xSym));
	double *B = REAL(GET_SLOT(val, Matrix_xSym));

	if(GPUFlag == 1) {
	
#ifdef HIPLAR_DBG	
		R_ShowMessage("DBG: Solving system of Ax = b, A = dpo, b = dge, using dpotrs_gpu;");
#endif
		double *d_A, *d_B;
		const char *uplo = uplo_P(Chol);
		magma_int_t NRHS = bdims[1];
		magma_int_t lda	 = adims[1];
		magma_int_t ldb  = bdims[0];
		magma_int_t N 	 = adims[0];
		cublasStatus retStatus;
		int i;

		/*if(uplo == "U")
			uplo = MagmaUpperStr;
		else if(uplo == "L")
			uplo = MagmaLowerStr;
		else		
			uplo = MagmaUpperStr;
		*/

		cublasAlloc(N * lda, sizeof(double), (void**)&d_A);
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation"));
		/********************************************/

		cublasAlloc(N * NRHS, sizeof(double), (void**)&d_B);	

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation"));
		/********************************************/

		cublasSetVector( N  * lda , sizeof(double), A, 1, d_A, 1);
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		cublasSetVector( ldb * NRHS, sizeof(double), B, 1, d_B, 1 );
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/

		magma_dpotrs_gpu(uplo[0], N, NRHS , d_A, lda, d_B, ldb, &info);

		cublasGetVector( ldb * NRHS, sizeof(double), d_B, 1, B, 1);
		
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer from Device"));
		/********************************************/
		
		cublasFree(d_A);
		cublasFree(d_B);
	}
	else {
	
#ifdef HIPLAR_DBG	
		R_ShowMessage("DBG: Solving system of Ax = b, A = dpo, b = dge, using dpotrs;");
#endif
		F77_CALL(dpotrs)(uplo_P(Chol), adims, bdims + 1, A , adims, B , bdims, &info);
	}
	if (info) {
		if(info > 0)
			error(_("the leading minor of order %d is not positive definite"),
					info);
		else /* should never happen! */
			error(_("Lapack routine %s returned error code %d"), "dpotrf", info);
	}
	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}

SEXP magma_dpoMatrix_matrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_MAGMA
    SEXP Chol = magma_dpoMatrix_chol(a),
	val = PROTECT(duplicate(b));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(getAttrib(b, R_DimSymbol)),
	info;

    if (!(isReal(b) && isMatrix(b)))
	error(_("Argument b must be a numeric matrix"));
    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    
		double *A = REAL(GET_SLOT(Chol, Matrix_xSym));
		double *B = REAL(val);
		const char *uplo = uplo_P(Chol);
		int N = bdims[1];
		//There is only a GPU interface for this call
		//so it will be the default setting if the GPU is on
		if(GPUFlag == 1) {
	
#ifdef HIPLAR_DBG	
			R_ShowMessage("DBG: Solving system of Ax = b, A = dpo, b = dge, using dpotrs_gpu;");
#endif
			double *d_A, *d_B;
			const char *uplo = uplo_P(Chol);
			magma_int_t NRHS = bdims[1];
			magma_int_t lda	 = adims[1];
			magma_int_t ldb  = bdims[0];
			magma_int_t N 	 = adims[0];
			cublasStatus retStatus;
			int i;
			cublasAlloc(N * lda, sizeof(double), (void**)&d_A);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory Allocation"));
			/********************************************/

			cublasAlloc(N * NRHS, sizeof(double), (void**)&d_B);	

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory Allocation"));
			/********************************************/

			cublasSetVector( N  * lda , sizeof(double), A, 1, d_A, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			cublasSetVector( ldb * NRHS, sizeof(double), B, 1, d_B, 1 );

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			magma_dpotrs_gpu(uplo[0], N, NRHS , d_A, lda, d_B, ldb, &info);

			cublasGetVector( ldb * NRHS, sizeof(double), d_B, 1, B, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer from Device"));
			/********************************************/

			cublasFree(d_A);
			cublasFree(d_B);
		}
		else {
		F77_CALL(dpotrs)(uplo_P(Chol), adims, bdims + 1,
		     REAL(GET_SLOT(Chol, Matrix_xSym)), adims,
		     REAL(val), bdims, &info);
		}
		// Error checking of MAGMA/LAPACK calls
		if (info) {
			if(info > 0)
				error(_("the leading minor of order %d is not positive definite"),
						info);
			else /* should never happen! */
				error(_("Lapack routine %s returned error code %d"), "dpotrf", info);
		}

		UNPROTECT(1);
    return val;
#endif
	return R_NilValue;
}
