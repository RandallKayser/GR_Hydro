#include "constants.h"

#ifndef INCLUDE_TIMESTEP
#define INCLUDE_TIMESTEP 1

void timestep_basic(double *Ustate_ptr, double *Pstate_ptr, double dt);

void RK2(double *Ustate_ptr, double *Pstate_ptr, double dt);
void RK3(double *Ustate_ptr, double *Pstate_ptr, double dt);
void RK3_TVD(double *Ustate_ptr, double *Pstate_ptr, double dt);
void RK4((double *Ustate_ptr, double *Pstate_ptr, double dt);

#endif