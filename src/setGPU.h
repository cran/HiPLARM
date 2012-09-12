#ifndef __GPU_FLAG__
#define __GPU_FLAG__

#include "Mutils.h"

extern int GPUFlag;
extern int Interface;

SEXP setGPU(SEXP, SEXP);
void interfaceSet(const char * );

#endif
