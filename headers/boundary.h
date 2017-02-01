#include "constants.h"

void enforce_bc(double *Ustate_ptr, double *Pstate_ptr, double x0, int x1type, int x2type, int x3type);
// 0 periodic
// 1 fixed_value
// 2 zero_gradient
// 3 bc_reflecting

void bc_periodic(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]);
void bc_fixed_value(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]);
void bc_zero_gradient(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]);
void bc_reflecting(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]);