#include <string.h>
#include "Mutils.h"
#include "hiplar_at.h"
#include "hiplar_dbg.h"
#include <R_ext/Print.h>
int hiplar_library = HIPLAR_USE_AUTO;


/* dgeMatrix */

int xover_dgeMatrix_norm                 = 0;
int xover_dgeMatrix_rcond                = 0;
int xover_dgeMatrix_crossprod            = 0;
int xover_dgeMatrix_dgeMatrix_crossprod  = 0;
int xover_dgeMatrix_matrix_crossprod     = 0;
int xover_dgeMatrix_LU                   = 0;
int xover_dgeMatrix_determinant          = 0;
int xover_dgeMatrix_solve                = 0;
int xover_dgeMatrix_matrix_solve         = 0;
int xover_dgeMatrix_matrix_mm            = 0;


/* dpoMatrix */

int xover_dpoMatrix_rcond                = 0;
int xover_dpoMatrix_solve                = 0;
int xover_dpoMatrix_matrix_solve         = 0;
int xover_dpoMatrix_dgeMatrix_solve      = 0;
int xover_dpoMatrix_chol                 = 0;


/* dsyMatrix */

int xover_dsyMatrix_matrix_mm            = 0;
int xover_dsyMatrix_norm                 = 0;


/* dtrMatrix */

int xover_dtrMatrix_solve                = 0;
int xover_dtrMatrix_chol2inv             = 0;
int xover_dtrMatrix_matrix_solve         = 0;
int xover_dtrMatrix_matrix_mm            = 0;
int xover_dtrMatrix_dtrMatrix_mm         = 0;



SEXP hiplarSet(SEXP var, SEXP val) {

	const char *pVar = CHAR(STRING_ELT(var,0));
    int tmpVal = asInteger(val);
		
	if(strcmp(pVar, "hiplar_library") == 0) {
        switch (tmpVal) {
            case HIPLAR_USE_PLASMA:
                hiplar_library = HIPLAR_USE_PLASMA;
                break;
            case HIPLAR_USE_MAGMA:
                hiplar_library = HIPLAR_USE_MAGMA;
                break;
            case HIPLAR_USE_AUTO:
                hiplar_library = HIPLAR_USE_AUTO;
                break;
            default:
                hiplar_library = HIPLAR_USE_AUTO;
        }

	} else if(strcmp(pVar, "xover_dgeMatrix_norm") == 0) {
        xover_dgeMatrix_norm = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_rcond") == 0) {
        xover_dgeMatrix_rcond = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_crossprod") == 0) {
        xover_dgeMatrix_crossprod = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_dgeMatrix_crossprod") == 0) {
        xover_dgeMatrix_dgeMatrix_crossprod = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_matrix_crossprod") == 0) {
        xover_dgeMatrix_matrix_crossprod = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_LU") == 0) {
        xover_dgeMatrix_LU = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_determinant") == 0) {
        xover_dgeMatrix_determinant = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_solve") == 0) {
        xover_dgeMatrix_solve = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_matrix_solve") == 0) {
        xover_dgeMatrix_matrix_solve = tmpVal;
	} else if(strcmp(pVar, "xover_dgeMatrix_matrix_mm") == 0) {
        xover_dgeMatrix_matrix_mm = tmpVal;

	} else if(strcmp(pVar, "xover_dpoMatrix_rcond") == 0) {
        xover_dpoMatrix_rcond = tmpVal;
	} else if(strcmp(pVar, "xover_dpoMatrix_solve") == 0) {
        xover_dpoMatrix_solve = tmpVal;
	} else if(strcmp(pVar, "xover_dpoMatrix_matrix_solve") == 0) {
        xover_dpoMatrix_matrix_solve = tmpVal;
	} else if(strcmp(pVar, "xover_dpoMatrix_dgeMatrix_solve") == 0) {
        xover_dpoMatrix_dgeMatrix_solve = tmpVal;
	} else if(strcmp(pVar, "xover_dpoMatrix_chol") == 0) {
        xover_dpoMatrix_chol = tmpVal;

	} else if(strcmp(pVar, "xover_dsyMatrix_matrix_mm") == 0) {
        xover_dsyMatrix_matrix_mm = tmpVal;
	} else if(strcmp(pVar, "xover_dsyMatrix_norm") == 0) {
        xover_dsyMatrix_norm = tmpVal;

	} else if(strcmp(pVar, "xover_dtrMatrix_solve") == 0) {
        xover_dtrMatrix_solve = tmpVal;
	} else if(strcmp(pVar, "xover_dtrMatrix_chol2inv") == 0) {
        xover_dtrMatrix_chol2inv = tmpVal;
	} else if(strcmp(pVar, "xover_dtrMatrix_matrix_solve") == 0) {
        xover_dtrMatrix_matrix_solve = tmpVal;
	} else if(strcmp(pVar, "xover_dtrMatrix_matrix_mm") == 0) {
        xover_dtrMatrix_matrix_mm = tmpVal;
	} else if(strcmp(pVar, "xover_dtrMatrix_dtrMatrix_mm") == 0) {
        xover_dtrMatrix_dtrMatrix_mm = tmpVal;
    } else {
        R_ShowMessage("ERROR: Unknow variable");
    }

    return R_NilValue;

}


SEXP hiplarShow() {

   Rprintf("HiPLAR library: ");
   switch (hiplar_library) {
        case HIPLAR_USE_PLASMA:
            Rprintf("PLASMA\n");
            break;
        case HIPLAR_USE_MAGMA:
            Rprintf("MAGMA\n");
            break;
        case HIPLAR_USE_AUTO:
            Rprintf("auto\n");
            break;
    }

    Rprintf("xover_dgeMatrix_norm                = %5d\n", xover_dgeMatrix_norm);
    Rprintf("xover_dgeMatrix_rcond               = %5d\n", xover_dgeMatrix_rcond);
    Rprintf("xover_dgeMatrix_crossprod           = %5d\n", xover_dgeMatrix_crossprod);
    Rprintf("xover_dgeMatrix_dgeMatrix_crossprod = %5d\n", xover_dgeMatrix_dgeMatrix_crossprod);
    Rprintf("xover_dgeMatrix_matrix_crossprod    = %5d\n", xover_dgeMatrix_matrix_crossprod);
    Rprintf("xover_dgeMatrix_LU                  = %5d\n", xover_dgeMatrix_LU);
    Rprintf("xover_dgeMatrix_determinant         = %5d\n", xover_dgeMatrix_determinant);
    Rprintf("xover_dgeMatrix_solve               = %5d\n", xover_dgeMatrix_solve);
    Rprintf("xover_dgeMatrix_matrix_solve        = %5d\n", xover_dgeMatrix_matrix_solve);
    Rprintf("xover_dgeMatrix_matrix_mm           = %5d\n", xover_dgeMatrix_matrix_mm);

    Rprintf("xover_dpoMatrix_rcond               = %5d\n", xover_dpoMatrix_rcond);
    Rprintf("xover_dpoMatrix_solve               = %5d\n", xover_dpoMatrix_solve);
    Rprintf("xover_dpoMatrix_matrix_solve        = %5d\n", xover_dpoMatrix_matrix_solve);
    Rprintf("xover_dpoMatrix_dgeMatrix_solve     = %5d\n", xover_dpoMatrix_dgeMatrix_solve);
    Rprintf("xover_dpoMatrix_chol                = %5d\n", xover_dpoMatrix_chol);

    Rprintf("xover_dsyMatrix_matrix_mm           = %5d\n", xover_dsyMatrix_matrix_mm);
    Rprintf("xover_dsyMatrix_norm                = %5d\n", xover_dsyMatrix_norm);

    Rprintf("xover_dtrMatrix_solve               = %5d\n", xover_dtrMatrix_solve);
    Rprintf("xover_dtrMatrix_chol2inv            = %5d\n", xover_dtrMatrix_chol2inv);
    Rprintf("xover_dtrMatrix_matrix_solve        = %5d\n", xover_dtrMatrix_matrix_solve);
    Rprintf("xover_dtrMatrix_matrix_mm           = %5d\n", xover_dtrMatrix_matrix_mm);
    Rprintf("xover_dtrMatrix_dtrMatrix_mm        = %5d\n", xover_dtrMatrix_dtrMatrix_mm);

    return R_NilValue;

}


SEXP hiplarWithPLASMA() {

#ifdef HIPLAR_WITH_PLASMA
    return ScalarLogical(TRUE);
#else
    return ScalarLogical(FALSE);
#endif

}


SEXP hiplarWithMAGMA() {

#ifdef HIPLAR_WITH_MAGMA
    return ScalarLogical(TRUE);
#else
    return ScalarLogical(FALSE);
#endif

}


