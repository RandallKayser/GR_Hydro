#include "constants.h"

void enforce_bc(double *Ustate_ptr, int x1type, int x2type, int x3type);
// 0 periodic
// 1 fixed_value
// 2 zero_gradient
// 3 bc_reflecting

void bc_periodic(int cell[3]);
void bc_fixed_value(int cell[3]);
void bc_zero_gradient(int cell[3]);
void bc_reflecting(int cell[3]);