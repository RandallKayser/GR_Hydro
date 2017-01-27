#include "constants.h"
#ifndef INCLUDE_INIT_FUNCS
#define INCLUDE_INIT_FUNCS 1

//primitive vars at x0 = 0
double rho_0(int cell[DIM_NUM-1]);
double p_0(int cell[DIM_NUM-1]);
double u0_0(int cell[DIM_NUM-1]);
double u1_0(int cell[DIM_NUM-1]);
double u2_0(int cell[DIM_NUM-1]);
double u3_0(int cell[DIM_NUM-1]);

void Pinit(double *Pstate_ptr);

double *U_ptr;
double *P_ptr;

#endif