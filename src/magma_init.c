#include <stdlib.h>
#include <stdio.h>
#include <R.h>

#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#include <cuda.h>
#include "minit.h"
CUdevice  dev;
CUcontext context;
#endif


int magma_initialise() {
#ifdef HIPLAR_WITH_MAGMA

    if(CUDA_SUCCESS != cuInit(0)) {
        R_ShowMessage("HiPLAR: cuInit failed\n" );
        exit(-1);
    }

    if(CUDA_SUCCESS != cuDeviceGet(&dev, 0)) {
        R_ShowMessage("HiPLAR: cuDeviceGet failed\n");
        exit(-1);
    }

    if(CUDA_SUCCESS != cuCtxCreate(&context, 0, dev)) {
        R_ShowMessage("HiPLAR: cuCtxCreate failed");
        exit(-1);
    }

    if(CUBLAS_STATUS_SUCCESS != cublasInit()) {
        R_ShowMessage("HiPLAR: cublasInit failed");
        exit(-1);
    }

    printout_devices();

#endif
	return 0;
}


int magma_deinit() {
#ifdef HIPLAR_WITH_MAGMA

    cuCtxDetach( context );
    cublasShutdown();

#endif
    return 0;
}
