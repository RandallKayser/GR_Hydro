#include "initfuncs.h"
#include "misc.h"

double rho_0(int cell[DIM_NUM-1]) {
	return pow(-1. / metric_full(0, 0, t_init, cell), -.5);
}


double p_0(int cell[DIM_NUM-1]) {
	return 10.0;
}


double u0_0(int cell[DIM_NUM-1]) {
	return 1.0;
}


double u1_0(int cell[DIM_NUM-1]) {
	return 0.0;
}


double u2_0(int cell[DIM_NUM-1]) {
	return 0.0;
}


double u3_0(int cell[DIM_NUM-1]) {
	return 0.0;
}


void Pinit(double *Pstate_ptr) {
	int cell[DIM_NUM-1];
	double *Pstate_local;
	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
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
