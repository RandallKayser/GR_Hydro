#include "../headers/constants.h"
#include "../headers/timestep.h"
#include "../headers/geometry.h"
#include "../headers/update.h"
#include "../headers/boundary.h"
#include "../headers/misc.h"
#include <stdlib.h>

void forward_euler(double *Ustate_ptr, double *U1state_ptr, double *Pstate_ptr, double *x0) {
	prim_to_cons(Ustate_ptr, Pstate_ptr, *x0);
	double dx;
	dx = calculate_dx_min(*x0);
	double dt = get_dt_max(Pstate_ptr, *x0, dx);

	enforce_bc(Ustate_ptr, Pstate_ptr, *x0, x1type, x2type, x3type);
	update(Ustate_ptr, Ustate_ptr, U1state_ptr, Pstate_ptr, 0., 1., dt, *x0);
	cons_to_prim(U1state_ptr, Pstate_ptr, *x0);
	
	*x0 += dt;
}


void RK2_TVD(double *Ustate_ptr, double *U1state_ptr, double *Pstate_ptr, double *x0) {
	prim_to_cons(Ustate_ptr, Pstate_ptr, *x0);
	double dx;
	dx = calculate_dx_min(*x0);
	double dt = get_dt_max(Pstate_ptr, *x0, dx);

	enforce_bc(Ustate_ptr, Pstate_ptr, *x0, x1type, x2type, x3type);
	update(Ustate_ptr, Ustate_ptr, U1state_ptr, Pstate_ptr, 0., 1., dt, *x0);
	cons_to_prim(U1state_ptr, Pstate_ptr, *x0);
	
	enforce_bc(U1state_ptr, Pstate_ptr, *x0, x1type, x2type, x3type);
	update(Ustate_ptr, U1state_ptr, Ustate_ptr, Pstate_ptr, .5, .5, dt, *x0);
	cons_to_prim(Ustate_ptr, Pstate_ptr, *x0);

	*x0 += dt;
}

void RK3_TVD(double *Ustate_ptr, double *U1state_ptr, double *Pstate_ptr, double *x0) {
	prim_to_cons(Ustate_ptr, Pstate_ptr, *x0);
	double dx;
	dx = calculate_dx_min(*x0);
	double dt = get_dt_max(Pstate_ptr, *x0, dx);
	
	enforce_bc(Ustate_ptr, Pstate_ptr, *x0, x1type, x2type, x3type);
	update(Ustate_ptr, Ustate_ptr, U1state_ptr, Pstate_ptr, 0., 1., dt, *x0);
	cons_to_prim(U1state_ptr, Pstate_ptr, *x0);
	
	enforce_bc(U1state_ptr, Pstate_ptr, *x0, x1type, x2type, x3type);
	update(Ustate_ptr, U1state_ptr, U1state_ptr, Pstate_ptr, .75, .25, dt, *x0);
	cons_to_prim(U1state_ptr, Pstate_ptr, *x0);
	
	enforce_bc(U1state_ptr, Pstate_ptr, *x0, x1type, x2type, x3type);
	update(Ustate_ptr, U1state_ptr, Ustate_ptr, Pstate_ptr, 1./3., 2./3., dt, *x0);
	cons_to_prim(Ustate_ptr, Pstate_ptr, *x0);

	*x0 += dt;
}

