#include "Mutils.h"
#include "Syms.h"

void minit() {

    Matrix_DimNamesSym = install("Dimnames");
    Matrix_DimSym = install("Dim");
    Matrix_diagSym = install("diag");
    Matrix_factorSym = install("factors");
    Matrix_iSym = install("i");
    Matrix_jSym = install("j");
    Matrix_lengthSym = install("length");
    Matrix_pSym = install("p");
    Matrix_permSym = install("perm");
    Matrix_uploSym = install("uplo");
    Matrix_xSym = install("x");

}
