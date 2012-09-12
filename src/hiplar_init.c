#include "minit.h"
#include "plasma_init.h"

void hiplar_init() {

	plasma_init();
    magma_initialise();
    minit();

}

void hiplar_deinit() {

	plasma_deinit();
    magma_deinit();

}
