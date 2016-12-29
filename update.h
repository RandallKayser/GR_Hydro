#include "constants.h"


double h(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);
double v_i(double *Pstate_ptr, int index, double x0, int cell[DIM_NUM-1]);
double W(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);
double D(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);
double E(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);
double S_k(double *Pstate_ptr, int index, double x0, int cell[DIM_NUM-1]);
void prim_to_cons(double *Ustate_ptr, double *Pstate_ptr, double x0);


double get_signals(int plus, double *Pstate_ptr, int dir, double x0, int cell[DIM_NUM-1]);
double get_dt_max(double *Pstate_ptr, double x0);
double* get_fluxes(double *Pstate_ptr, double* (*riemann_input)(double*, double, int[DIM_NUM-1]));
void update(double *Ustate_ptr, double *Pstate_ptr, double x0);

double* piecewise_constant(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);
double* plm(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);
double* ppm(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]);

void cons_to_prim(double *Ustate_ptr, double *Pstate_ptr, double x0);








