#include <string.h>
#include "setGPU.h"
#include <R_ext/Print.h>
int GPUFlag = 1;
int Interface = 0;

SEXP setGPU(SEXP x, SEXP interface) {

	int *tmp = INTEGER(x);
	const char *val = CHAR(STRING_ELT(interface,0));
		
	if(*tmp > 0) {
			GPUFlag = 1;
			Rprintf("GPU Active in %s mode\n", val);
		}
		else {
			GPUFlag = 0;
			Rprintf("GPU Inactive\n");
		}
	if(val) {
		interfaceSet(val);
	} else
		Interface = 0;

	return x;
}

void interfaceSet(const char *val) {


	if(GPUFlag == 1) {
		if(strcmp(val,"GPU") == 0) 
			Interface = 1;
		else if(strcmp(val,"CPU") == 0)
			Interface = 0;
		else 
			Interface = 0;
	}

}



