#include "constants.h"

int bc_type(int dir);
// 0 dirichlet
// 1 von neumann
// 2 periodic
void bc_periodic();
void bc_dirichlet();
void bc_freeflow();