#include <math.h>
#include "initfuncs.h"

double rho_0(int cell[DIM_NUM-1]) {
	return 1.0;
}


double p_0(int cell[DIM_NUM-1]) {
	return 1.0;
}


double u0_0(int cell[DIM_NUM-1]) {
	return 1.0;
}


double u1_0(int cell[DIM_NUM-1]) {
	return 1.0;
}


double u2_0(int cell[DIM_NUM-1]) {
	return 1.0;
}


double u3_0(int cell[DIM_NUM-1]) {
	return 1.0;
}


void Pinit(double *Pstate_ptr) {
	int cell[DIM_NUM-1];

	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;

				*(Pstate_ptr + P_offset(cell, 0)) = rho_0(cell);
				*(Pstate_ptr + P_offset(cell, 1)) = p_0(cell);
				*(Pstate_ptr + P_offset(cell, 2)) = u0_0(cell);
				*(Pstate_ptr + P_offset(cell, 3)) = u1_0(cell);
				*(Pstate_ptr + P_offset(cell, 4)) = u2_0(cell);		
				*(Pstate_ptr + P_offset(cell, 5)) = u3_0(cell);
			}
		}
	}
}
