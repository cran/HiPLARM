#include <stdlib.h>
#include <string.h>
#include <R_ext/Print.h>
#ifdef HIPLAR_WITH_PLASMA
#include <plasma.h>
#endif


int R_PLASMA_MAJOR;
int R_PLASMA_MINOR;
int R_PLASMA_MICRO;

int R_PLASMA_NUM_THREADS;
int R_PLASMA_SCHED;


int plasma_init() {
#ifdef HIPLAR_WITH_PLASMA
    int info;

	R_PLASMA_NUM_THREADS = 1;

	if (getenv("R_PLASMA_NUM_THREADS") != NULL) {
		R_PLASMA_NUM_THREADS = atoi(getenv("R_PLASMA_NUM_THREADS"));
	} else {
		  Rprintf("The envirnment variable R_PLASMA_NUM_THREADS is not set.\n");
		  Rprintf("Using one thread for PLASMA.\nPlease set R_PLASMA_NUM_THREADS to the number of actual cores.\n\n");
			//printf("\nThe envirnment variable R_PLASMA_NUM_THREADS is not set.\n");
			//printf("Using one thread for PLASMA.\nPlease set R_PLASMA_NUM_THREADS to the number of actual cores.\n\n");v
    }

    /* Init PLASMA */
    info = PLASMA_Init(R_PLASMA_NUM_THREADS);

    if ((getenv("R_PLASMA_SCHED") != NULL) &&
        (strcmp(getenv("R_PLASMA_SCHED"), "STATIC") == 0)) {
        PLASMA_Set( PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );
        R_PLASMA_SCHED = 0;
    } else {
        PLASMA_Set( PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
        R_PLASMA_SCHED = 1;
    }

    PLASMA_Version(&R_PLASMA_MAJOR, &R_PLASMA_MINOR, &R_PLASMA_MICRO);

    if ((PLASMA_VERSION_MAJOR != R_PLASMA_MAJOR) ||
        (PLASMA_VERSION_MINOR != R_PLASMA_MINOR) ||
        (PLASMA_VERSION_MICRO != R_PLASMA_MICRO)) {
        Rprintf("ERROR: PLASMA version mismatch\n");
	  Rprintf("ERROR: PLASMA version mismatch %d.%d.%d %d.%d.%d\n",
    PLASMA_VERSION_MAJOR, PLASMA_VERSION_MINOR, PLASMA_VERSION_MICRO,
    R_PLASMA_MAJOR, R_PLASMA_MINOR, R_PLASMA_MICRO);
    }

	return(info);
#endif
	return 0;
}


int plasma_deinit() {
#ifdef HIPLAR_WITH_PLASMA

    /* Deinit PLASMA */
    return PLASMA_Finalize();

#endif
    return 0;
}
