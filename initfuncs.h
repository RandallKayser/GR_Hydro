#include "constants.h"

//primitive vars at x0 = 0
double rho_0(int cell[DIM_NUM-1]);
double p_0(int cell[DIM_NUM-1]);
double u0_0(int cell[DIM_NUM-1]);
double u1_0(int cell[DIM_NUM-1]);
double u2_0(int cell[DIM_NUM-1]);
double u3_0(int cell[DIM_NUM-1]);

void Pinit(double *Pstate_ptr);

double U[x1cellnum * x2cellnum * x3cellnum * (1+DIM_NUM)];
double Prims[x1cellnum * x2cellnum * x3cellnum * (2+DIM_NUM)];
double *U_ptr;
double *Prims_ptr;

U_ptr = U;
Prims_ptr = Prims;