#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#endif
#include "magma_dgeMatrix.h"
#include "setGPU.h"
#include "hiplar_dbg.h"

static
double magma_get_norm(SEXP obj, const char *typstr)
{
#ifdef HIPLAR_WITH_MAGMA
	if(any_NA_in_x(obj))
		return NA_REAL;
	else {
		char typnm[] = {'\0', '\0'};
		int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
		double *work = (double *) NULL;

		typnm[0] = La_norm_type(typstr);
		if (*typnm == 'I') {
			work = (double *) R_alloc(dims[0], sizeof(double));
			if(GPUFlag == 1 && (dims[0] % 64 == 0) && (dims[1] % 64 == 0)) {
#ifdef HIPLAR_DBG
	R_ShowMessage("DBG: Getting norm using magmablas_dlange");
#endif
				double *d_work, *d_A, *A, val;
				A = REAL(GET_SLOT(obj, Matrix_xSym));
				cublasAlloc(dims[0] * dims[1], sizeof(double), (void**)&d_A);
				cublasAlloc(dims[0], sizeof(double), (void**)&d_work);
				cublasSetVector(dims[0] * dims[1], sizeof(double), A, 1, d_A, 1);
				val = magmablas_dlange(*typstr, dims[0], dims[1], d_A, dims[0], d_work);
				cudaFree(d_A);
				cudaFree(d_work);
				return val;
			}

		}
		return F77_CALL(dlange)(typstr, dims, dims+1,
				REAL(GET_SLOT(obj, Matrix_xSym)),
				dims, work);
	}
#endif

	return 0.0;
}

SEXP magma_dgeMatrix_norm(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_MAGMA
	    return ScalarReal(magma_get_norm(obj, CHAR(asChar(type))));
#endif
    return R_NilValue;
}


SEXP magma_dgeMatrix_LU_(SEXP x, Rboolean warn_sing)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = get_factors(x, "LU");
	int *dims, npiv, info;

	if (val != R_NilValue) {
//		R_ShowMessage("already in slot");	/* nothing to do if it's there in 'factors' slot */
		return val;
	}

	dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
	if (dims[0] < 1 || dims[1] < 1)
		error(_("Cannot factor a matrix with zero extents"));
	npiv = (dims[0] < dims[1]) ? dims[0] : dims[1];
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("denseLU")));
	slot_dup(val, x, Matrix_xSym);
	slot_dup(val, x, Matrix_DimSym);
	double *h_R = REAL(GET_SLOT(val, Matrix_xSym));
	int *ipiv = INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, npiv));
	
	if(GPUFlag == 0){
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: LU decomposition using dgetrf;");
#endif
		F77_CALL(dgetrf)(dims, dims + 1, h_R,
				dims,
				ipiv,
				&info);
	}
	else if(GPUFlag == 1 && Interface == 0){
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: LU decomposition using magma_dgetrf;");
#endif
			magma_dgetrf(dims[0], dims[1], h_R, dims[0], ipiv, &info);
		}
		else if(GPUFlag == 1 && Interface == 1) {
			
#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: LU decomposition using magma_dgetrf_gpu;");	
#endif
			double *d_A;
			int N2 = dims[0] * dims[1];
			cublasStatus retStatus;

			cublasAlloc( N2 , sizeof(double), (void**)&d_A);
			
			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory Allocation"));
			/********************************************/

			cublasSetVector(N2, sizeof(double), h_R, 1, d_A, 1);
			
			/* Error Checking */
			retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Date Transfer to Device"));
			/********************************************/

			magma_dgetrf_gpu(dims[0],dims[1], d_A, dims[0], ipiv,  &info);
			
			cublasGetVector( N2, sizeof(double), d_A, 1, h_R, 1);		
			
			/* Error Checking */
			retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Date Transfer from Device"));
			/********************************************/
			
				cublasFree(d_A);
		
			/* Error Checking */
			retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error freeing data"));
			/********************************************/
		}
		else
			error(_("MAGMA/LAPACK/Interface Flag not defined correctly"));
		
	if (info < 0)
		error(_("Lapack routine %s returned error code %d"), "dgetrf", info);
	else if (info > 0 && warn_sing)
		warning(_("Exact singularity detected during LU decomposition: %s, i=%d."),
				"U[i,i]=0", info);
	UNPROTECT(1);

	return set_factors(x, val, "LU");
#endif

	    return R_NilValue;
}

SEXP magma_dgeMatrix_LU(SEXP x, SEXP warn_singularity)
{
#ifdef HIPLAR_WITH_MAGMA
	    return magma_dgeMatrix_LU_(x, asLogical(warn_singularity));
#endif	
	    return R_NilValue;
}

/* Solvers for dgeMatrices using LUdecomp */

SEXP magma_dgeMatrix_solve(SEXP a)
{
#ifdef HIPLAR_WITH_MAGMA
    /*  compute the 1-norm of the matrix, which is needed
	later for the computation of the reciprocal condition number. */
    double aNorm = magma_get_norm(a, "1");

    /* the LU decomposition : */
		/* Given that we may be performing this operation
		 * on the GPU we may put in an optimisation here
		 * where if we call the LU solver we, we do not require
		 * the decomposition to be transferred back to CPU. This is TODO
		 */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	lu = magma_dgeMatrix_LU_(a, TRUE);
    int *dims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*pivot = INTEGER(GET_SLOT(lu, Matrix_permSym));

    /* prepare variables for the dgetri calls */
    double *x, tmp;
    int	info, lwork = -1;

    if (dims[0] != dims[1]) error(_("Solve requires a square matrix"));
    slot_dup(val, lu, Matrix_xSym);
    x = REAL(GET_SLOT(val, Matrix_xSym));
    slot_dup(val, lu, Matrix_DimSym);
		int N2 = dims[0] * dims[0];

    if(dims[0]) /* the dimension is not zero */
    {
			/* is the matrix is *computationally* singular ? */
			double rcond;
			F77_CALL(dgecon)("1", dims, x, dims, &aNorm, &rcond,
					(double *) R_alloc(4*dims[0], sizeof(double)),
					(int *) R_alloc(dims[0], sizeof(int)), &info);
			if (info)
				error(_("error [%d] from Lapack 'dgecon()'"), info);
			if(rcond < DOUBLE_EPS)
				error(_("Lapack dgecon(): system computationally singular, reciprocal condition number = %g"),
						rcond);

			/* only now try the inversion and check if the matrix is *exactly* singular: */
			// This is also a work space query. This is not an option in magma

			F77_CALL(dgetri)(dims, x, dims, pivot, &tmp, &lwork, &info);
			lwork = (int) tmp;
			
			if( GPUFlag == 0){
				

				F77_CALL(dgetri)(dims, x, dims, pivot,
						(double *) R_alloc((size_t) lwork, sizeof(double)),
						&lwork, &info);

#ifdef HIPLAR_DBG
				R_ShowMessage("DBG: Solve using LU using dgetri;");
#endif
			}
			else if(GPUFlag == 1) {
				
				double *d_x, *dwork, *dtmp; 
				cublasStatus retStatus;			
	
#ifdef HIPLAR_DBG
				R_ShowMessage("Solve using LU using magma_dgetri;");
#endif
				cublasAlloc(N2, sizeof(double), (void**)&d_x);

				//cublasAlloc(N2 , sizeof(double), (void**)&dtmp);
				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Memory Allocation on Device"));
				/********************************************/

				cublasSetVector( N2, sizeof(double), x, 1, d_x, 1);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data Transfer to Device"));
				/********************************************/
				lwork = dims[0] * magma_get_dgetri_nb( dims[0] );
				
					cublasAlloc(lwork, sizeof(double), (void**)&dwork);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Memory Allocation on Device"));
				/********************************************/

				magma_dgetri_gpu(dims[0], d_x, dims[0], pivot, dwork , lwork, &info);

				cublasGetVector(N2, sizeof(double), d_x, 1, x, 1);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data From to Device"));
				/********************************************/

				cublasFree(dwork);
				cublasFree(d_x);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error freeing memory"));
				/********************************************/
					
			}
			else
				error(_("GPUFlag not set correctly"));

			if (info)
				error(_("Lapack routine dgetri: system is exactly singular"));
    }
    UNPROTECT(1);
    return val;
#endif
	    return R_NilValue;
}

SEXP magma_dgeMatrix_matrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b)),
			 lu = PROTECT(magma_dgeMatrix_LU_(a, TRUE));
	int *adims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
			*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
	int info, n = bdims[0], nrhs = bdims[1];



	if (*adims != *bdims || bdims[1] < 1 || *adims < 1 || *adims != adims[1])
		error(_("Dimensions of system to be solved are inconsistent"));

	double *A = REAL(GET_SLOT(lu, Matrix_xSym));
	double *B  = REAL(GET_SLOT(val, Matrix_xSym));
	int *ipiv = INTEGER(GET_SLOT(lu, Matrix_permSym));

	if(GPUFlag == 0) {
		F77_CALL(dgetrs)("N", &n, &nrhs, A, &n, ipiv, B, &n, &info);	
	
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Solve using LU using dgetrs;");
#endif
	}else if(GPUFlag == 1) {
		
	
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Solve using LU using magma_dgetrs;");
#endif		
		double *d_A, *d_B;
		cublasStatus retStatus;

		cublasAlloc(adims[0] * adims[1], sizeof(double), (void**)&d_A);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation of A on Device"));
		/********************************************/


		cublasAlloc(n * nrhs, sizeof(double), (void**)&d_B);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation of b on Device"));
		/********************************************/



		cublasSetVector(adims[0] * adims[1], sizeof(double), A, 1, d_A, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Transferring data to advice"));
		/********************************************/

		cublasSetVector(n * nrhs, sizeof(double), B, 1, d_B, 1);

		magma_dgetrs_gpu( 'N', n, nrhs, d_A, n, ipiv, d_B, n, &info );

		cublasGetVector(n * nrhs, sizeof(double), d_B, 1, B, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Transferring from to advice"));
		/********************************************/

		cublasFree(d_A);
		cublasFree(d_B);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in freeing data"));
		/********************************************/

		
	}
	if (info)
		error(_("Lapack routine dgetrs: system is exactly singular"));
	UNPROTECT(2);
	return val;
#endif
	    return R_NilValue;
}


SEXP magma_dgeMatrix_matrix_mm(SEXP a, SEXP bP, SEXP right)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP b = PROTECT(mMatrix_as_dgeMatrix(bP)),
	     val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
	int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	    *bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	    *cdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
	double one = 1.0, zero = 0.0;

	if (asLogical(right)) {
		int m = bdims[0], n = adims[1], k = bdims[1];
		if (adims[0] != k)
			error(_("Matrices are not conformable for multiplication"));
		cdims[0] = m; cdims[1] = n;
		if (m < 1 || n < 1 || k < 1) {
			// 		This was commented out
				    error(_("Matrices with zero extents cannot be multiplied")); 
			ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n);
		} else {
			double *B = REAL(GET_SLOT(b, Matrix_xSym));
			double *A = REAL(GET_SLOT(a, Matrix_xSym));
			double *C = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n));
			//TODO add magma here too
			if(GPUFlag == 1) {
				double *d_A, *d_B, *d_C;
				cublasStatus retStatus;


#ifdef HIPLAR_DBG
				R_ShowMessage("DBG: Performing matrix multiplication with Right = true using magmablas_dgemm");
#endif
				cublasAlloc(n * k, sizeof(double), (void**)&d_A);
			
				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Memory Allocation"));
				/********************************************/

				cublasAlloc(m * k, sizeof(double), (void**)&d_B);	

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

				cublasSetVector( n  * k , sizeof(double), A, 1, d_A, 1);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data Transfer to Device"));
				/********************************************/

				cublasSetVector( m * k, sizeof(double), B, 1, d_B, 1 );

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data Transfer to Device"));
				/********************************************/

				// ******** magmablas_dgemm call Here **
				magmablas_dgemm('N', 'N', m, n, k, one, d_B, m, d_A, k, zero, d_C,  m);

				/* Error Checking */
			//	retStatus = cublasGetError ();
			//	if (retStatus != CUBLAS_STATUS_SUCCESS) 
			//		error(_("CUBLAS: Error in routine"));
				/********************************************/


				cublasGetVector( m * n , sizeof(double), d_C, 1, C, 1);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data Transfer from Device"));
				/********************************************/

				cublasFree(d_A);
				cublasFree(d_B);
				cublasFree(d_C);
			}
			else {
	
#ifdef HIPLAR_DBG
				R_ShowMessage("DBG: Performing matrix multiplication using dgemm with right = TRUE");
#endif
				F77_CALL(dgemm) ("N", "N", &m, &n, &k, &one,
						B, &m, A , &k, &zero, C , &m);
			}
		}
	} else {
		int m = adims[0], n = bdims[1], k = adims[1];
		double *A = REAL(GET_SLOT(a, Matrix_xSym));
		double *B = REAL(GET_SLOT(b, Matrix_xSym));
		double *C = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n));

		if (bdims[0] != k)
			error(_("Matrices are not conformable for multiplication"));
		cdims[0] = m; cdims[1] = n;
		if (m < 1 || n < 1 || k < 1) {
			//		This was commented out
			error(_("Matrices with zero extents cannot be multiplied")); 
			ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n);
		} else {
			if(GPUFlag == 1) {

				double *d_A, *d_B, *d_C;
				cublasStatus retStatus;

	
#ifdef HIPLAR_DBG
				R_ShowMessage("DBG: Performing matrix multiplication using magmablas_dgemm");
#endif			
				cublasAlloc(m * k, sizeof(double), (void**)&d_A);
			
				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Memory Allocation"));
				/********************************************/

				cublasAlloc(n * k, sizeof(double), (void**)&d_B);	

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

				cublasSetVector( m  * k , sizeof(double), A, 1, d_A, 1);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data Transfer to Device"));
				/********************************************/

				cublasSetVector( n * k, sizeof(double), B, 1, d_B, 1 );

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data Transfer to Device"));
				/********************************************/

				// ******** magmablas_dgemm call Here **
				magmablas_dgemm('N', 'N', m, n, k, one, d_A, m, d_B, k, zero, d_C,  m);

				cublasGetVector( m * n , sizeof(double), d_C, 1, C, 1);

				/* Error Checking */
				retStatus = cublasGetError ();
				if (retStatus != CUBLAS_STATUS_SUCCESS) 
					error(_("CUBLAS: Error in Data Transfer from Device"));
				/********************************************/

				cublasFree(d_A);
				cublasFree(d_B);
				cublasFree(d_C);
				
			}
			else {
	
#ifdef HIPLAR_DBG
				R_ShowMessage("DBG: Performing matrix multiplication using dgemm");
#endif
				F77_CALL(dgemm) ("N", "N", &m, &n, &k, &one,
						A, &m,
						B, &k, &zero,
						C,
						&m);	

			}
		}
	}
	ALLOC_SLOT(val, Matrix_DimNamesSym, VECSXP, 2);
	UNPROTECT(2);
	return val;
#endif
	return R_NilValue;
}

SEXP magma_dgeMatrix_crossprod(SEXP x, SEXP trans)
{
#ifdef HIPLAR_WITH_MAGMA
	int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x) */
	SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix"))),
	     nms = VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym), tr ? 0 : 1),
	     vDnms = ALLOC_SLOT(val, Matrix_DimNamesSym, VECSXP, 2);
	int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	    *vDims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
	int k = tr ? Dims[1] : Dims[0], n = tr ? Dims[0] : Dims[1];
	double *vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n * n)),
		 one = 1.0, zero = 0.0;
	double *A =  REAL(GET_SLOT(x, Matrix_xSym));
	AZERO(vx, n * n);
	SET_SLOT(val, Matrix_uploSym, mkString("U"));
	ALLOC_SLOT(val, Matrix_factorSym, VECSXP, 0);
	vDims[0] = vDims[1] = n;
	SET_VECTOR_ELT(vDnms, 0, duplicate(nms));
	SET_VECTOR_ELT(vDnms, 1, duplicate(nms));
	if(n && GPUFlag == 1) {

#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Performing crossproduct using cublasDsyrk");
#endif
		cublasStatus retStatus;
		cublasHandle_t handle;
		double *d_A, *d_C;

		/*retStatus = cublasCreate(&handle);
		  if ( retStatus != CUBLAS_STATUS_SUCCESS )		
		  error(_("CUBLAS initialisation failed"));
		  */

		cublasAlloc(n * k, sizeof(double), (void**)&d_A);
		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation"));
		/********************************************/

		cublasAlloc(n * n, sizeof(double), (void**)&d_C);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Memory Allocation"));
		/********************************************/

		cublasSetVector( n  * k , sizeof(double), A, 1, d_A, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/
		
		//cublasSetVector( n  * n , sizeof(double), vx, 1, d_C, 1);
		
		/* Error Checking */
		//retStatus = cublasGetError ();
		//if (retStatus != CUBLAS_STATUS_SUCCESS) 
		//	error(_("CUBLAS: Error in Data Transfer to Device"));
		/********************************************/


		cublasDsyrk('U' , tr ? 'N' : 'T', n, k, one, d_A, Dims[0], zero, d_C, n);

		cublasGetVector( n * n , sizeof(double), d_C, 1, vx, 1);

		/* Error Checking */
		retStatus = cublasGetError ();
		if (retStatus != CUBLAS_STATUS_SUCCESS) 
			error(_("CUBLAS: Error in Data Transfer from Device"));
		/********************************************/

		cublasFree(d_A);
		cublasFree(d_C);

	} else if(n){
	
#ifdef HIPLAR_DBG
		R_ShowMessage("DBG: Performing cross prod with dsyrk");
#endif
		F77_CALL(dsyrk)("U", tr ? "N" : "T", &n, &k, &one, A, Dims,
				&zero, vx, &n);
	}

	SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}

SEXP magma_dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
#ifdef HIPLAR_WITH_MAGMA
	int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */
	SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
	int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	    *yDims = INTEGER(GET_SLOT(y, Matrix_DimSym)),
	    *vDims;
	int m  = xDims[!tr],  n = yDims[!tr];/* -> result dim */
	int xd = xDims[ tr], yd = yDims[ tr];/* the conformable dims */
	double one = 1.0, zero = 0.0;
	double *A = REAL(GET_SLOT(x, Matrix_xSym));
	double *B = REAL(GET_SLOT(y, Matrix_xSym));

	SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
	SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
	vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
	if (xd > 0 && yd > 0 && n > 0 && m > 0) {
		if (xd != yd)
			error(_("Dimensions of x and y are not compatible for %s"),
					tr ? "tcrossprod" : "crossprod");
		vDims[0] = m; vDims[1] = n;
		SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
		double *C = REAL(GET_SLOT(val, Matrix_xSym));

		if(GPUFlag == 1){

			double *d_A, *d_B, *d_C;
			cublasStatus retStatus;
			int i;

#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Performing dge/dge crossprod using magmablas_dgemm");
#endif
			cublasAlloc(m * xd, sizeof(double), (void**)&d_A);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory Allocation"));
			/********************************************/

			cublasAlloc(n * xd, sizeof(double), (void**)&d_B);	

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

			cublasSetVector( m  * xd , sizeof(double), A, 1, d_A, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			cublasSetVector( xd * n, sizeof(double), B, 1, d_B, 1 );

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			//cublasSetVector( m * n, sizeof(double), C, 1, d_C, 1 );

			/* Error Checking */
			//retStatus = cublasGetError ();
			//if (retStatus != CUBLAS_STATUS_SUCCESS) 
			//	error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/


			// ******** magmablas_dgemm call Here **
			magmablas_dgemm( tr ? 'N' : 'T', tr ? 'T' : 'N', m, n, xd, one, d_A, xDims[0], d_B, yDims[0], zero, d_C,  m);

			cublasGetVector( m * n , sizeof(double), d_C, 1, C, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer from Device"));
			/********************************************/

			cublasFree(d_A);
			cublasFree(d_B);
			cublasFree(d_C);
		
		}
		else{
	
#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Performing dge/dge cross prod with dgemm");
#endif
			F77_CALL(dgemm)(tr ? "N" : "T", tr ? "T" : "N", &m, &n, &xd, &one, A , xDims, B , yDims, &zero, C, &m);
		}
	}
	UNPROTECT(1);
	return val;
#endif
	return R_NilValue;
}

SEXP magma_dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
#ifdef HIPLAR_WITH_MAGMA
	int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */
	SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
	int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	    *yDims = INTEGER(getAttrib(y, R_DimSymbol)),
	    *vDims, nprot = 1;
	int m  = xDims[!tr],  n = yDims[!tr];/* -> result dim */
	int xd = xDims[ tr], yd = yDims[ tr];/* the conformable dims */
	double one = 1.0, zero = 0.0;

	if (isInteger(y)) {
		y = PROTECT(coerceVector(y, REALSXP));
		nprot++;
	}
	if (!(isMatrix(y) && isReal(y)))
		error(_("Argument y must be a numeric matrix"));
	SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
	SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
	vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
	if (xd > 0 && yd > 0 && n > 0 && m > 0) {
		if (xd != yd)
			error(_("Dimensions of x and y are not compatible for %s"),
					tr ? "tcrossprod" : "crossprod");
		vDims[0] = m; vDims[1] = n;
		SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));

		double *A = REAL(GET_SLOT(x, Matrix_xSym));
		double *B = REAL(y);
		double *C = REAL(GET_SLOT(val, Matrix_xSym));

		if(GPUFlag == 1) {
			
			double *d_A, *d_B, *d_C;
			cublasStatus retStatus;
			int i;

#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Performing dge/matrix crossprod using magmablas_dgemm");
#endif
			cublasAlloc(m * xd, sizeof(double), (void**)&d_A);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Memory Allocation"));
			/********************************************/

			cublasAlloc(n * xd, sizeof(double), (void**)&d_B);	

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

			cublasSetVector( m  * xd , sizeof(double), A, 1, d_A, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			cublasSetVector( xd * n, sizeof(double), B, 1, d_B, 1 );

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/

			cublasSetVector( m * n, sizeof(double), C, 1, d_C, 1 );

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer to Device"));
			/********************************************/


			// ******** magmablas_dgemm call Here **
			magmablas_dgemm( tr ? 'N' : 'T', tr ? 'T' : 'N', m, n, xd, one, d_A, xDims[0], d_B, yDims[0], zero, d_C,  m);

			cublasGetVector( m * n , sizeof(double), d_C, 1, C, 1);

			/* Error Checking */
			retStatus = cublasGetError ();
			if (retStatus != CUBLAS_STATUS_SUCCESS) 
				error(_("CUBLAS: Error in Data Transfer from Device"));
			/********************************************/

			cublasFree(d_A);
			cublasFree(d_B);
			cublasFree(d_C);
			
		}
		else {

#ifdef HIPLAR_DBG
			R_ShowMessage("DBG: Performing dge/matrix cross prod with dgemm");
#endif
			F77_CALL(dgemm)(tr ? "N" : "T", tr ? "T" : "N", &m, &n, &xd, &one,
					A , xDims,
					B , yDims,
					&zero, C, &m);
		}
	}
	UNPROTECT(nprot);
	return val;
#endif
	return R_NilValue;
}

SEXP magma_dgeMatrix_determinant(SEXP x, SEXP logarithm)
{
#ifdef HIPLAR_WITH_MAGMA

#ifdef HIPLAR_DBG
	R_ShowMessage("DBG: Entering magma_dgeMatrix_determinant");
#endif

	int lg = asLogical(logarithm);
	int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	    n = dims[0], sign = 1;
	double modulus = lg ? 0. : 1; /* initialize; = result for n == 0 */

	if (n != dims[1])
		error(_("Determinant requires a square matrix"));
	if (n > 0) {
		SEXP lu = magma_dgeMatrix_LU_(x, /* do not warn about singular LU: */ FALSE);
		int i, *jpvt = INTEGER(GET_SLOT(lu, Matrix_permSym));
		double *luvals = REAL(GET_SLOT(lu, Matrix_xSym));

		for (i = 0; i < n; i++) if (jpvt[i] != (i + 1)) sign = -sign;
		if (lg) {
			for (i = 0; i < n; i++) {
				double dii = luvals[i*(n + 1)]; /* ith diagonal element */
				modulus += log(dii < 0 ? -dii : dii);
				if (dii < 0) sign = -sign;
			}
		} else {
			for (i = 0; i < n; i++)
				modulus *= luvals[i*(n + 1)];
			if (modulus < 0) {
				modulus = -modulus;
				sign = -sign;
			}
		}
	}
	return as_det_obj(modulus, lg, sign);
#endif
	return R_NilValue;
}

SEXP magma_dgeMatrix_rcond(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_MAGMA
	SEXP LU = PROTECT(magma_dgeMatrix_LU_(obj, FALSE));/* <- not warning about singularity */
	char typnm[] = {'\0', '\0'};
	int *dims = INTEGER(GET_SLOT(LU, Matrix_DimSym)), info;
	double anorm, rcond;

	if (dims[0] != dims[1] || dims[0] < 1) {
		UNPROTECT(1);
		error(_("rcond requires a square, non-empty matrix"));
	}
	typnm[0] = La_rcond_type(CHAR(asChar(type)));
	anorm = magma_get_norm(obj, typnm);
	F77_CALL(dgecon)(typnm,
			dims, REAL(GET_SLOT(LU, Matrix_xSym)),
			dims, &anorm, &rcond,
			(double *) R_alloc(4*dims[0], sizeof(double)),
			(int *) R_alloc(dims[0], sizeof(int)), &info);
	UNPROTECT(1);
	return ScalarReal(rcond);
#endif
	return R_NilValue;
}
