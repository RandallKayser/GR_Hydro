#include "headers/initfuncs.h"
#include "headers/misc.h"
#include "headers/geometry.h"
#include <math.h>

double rho_0(int cell[DIM_NUM-1]) {
	return 1.;
}


double p_0(int cell[DIM_NUM-1]) {
	return 10.;
}


double u0_0(int cell[DIM_NUM-1]) {
	return pow(-metric_full(0, 0, t_init, cell), -.5);
}


double u1_0(int cell[DIM_NUM-1]) {
	return 0.;
}


double u2_0(int cell[DIM_NUM-1]) {
	return 0.;
}


double u3_0(int cell[DIM_NUM-1]) {
	return 0.;
}


void Pinit(double *Pstate_ptr) {
	int cell[DIM_NUM-1];
	double *Pstate_local;
	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;

				Pstate_local = Pstate_ptr + P_offset(cell, 0);

				Pstate_local[0] = rho_0(cell);
				Pstate_local[1] = p_0(cell);
				Pstate_local[2] = u0_0(cell);
				Pstate_local[3] = u1_0(cell);
				Pstate_local[4] = u2_0(cell);		
				Pstate_local[5] = u3_0(cell);
			}
		}
	}
}
