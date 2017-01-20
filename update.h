#include "constants.h"
#include "timestep.h"

double h(double *Pstate_local, double x0, int cell[DIM_NUM-1]);
double v_i(double *Pstate_local, int index, double x0, int cell[DIM_NUM-1]);
double W(double *Pstate_local, double x0, int cell[DIM_NUM-1]);
double D(double *Pstate_local, double x0, int cell[DIM_NUM-1]);
double E(double *Pstate_local, double x0, int cell[DIM_NUM-1]);
double S_k(double *Pstate_local, int index, double x0, int cell[DIM_NUM-1]);
double T_mu_nu(double* Pstate_local, int dir1, int dir2, double x0, int cell[DIM_NUM-1]);
void prim_to_cons(double *Ustate_ptr, double *Pstate_ptr, double x0);


double get_signals(int plus, double *Pstate_ptr, int dir, double x0, int cell[DIM_NUM-1]);
double get_dt_max(double *Pstate_ptr, double x0);
double* fluxes(double *Pstate_ptr, double* (*flux)(double*, double* (*)(double*, double int[DIM_NUM-1])));
double* sources(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);
double* rhlle(double *Pstate_ptr, int dir, double x0, int cell[DIM_NUM-1]);

void update(double *Ustate_ptr, double *Pstate_ptr, double x0);

double* piecewise_constant(double *Pstate_ptr, int right_bias, int dir, double x0, int cell[DIM_NUM-1]);
double* plm(double *Pstate_ptr, int right_bias, int dir, double x0, int cell[DIM_NUM-1]);
double* ppm(double *Pstate_ptr,  int right_bias, int dir, double x0, int cell[DIM_NUM-1]);

void cons_to_prim(double *Ustate_ptr, double *Pstate_ptr, double x0);








