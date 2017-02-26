#include "../../headers/constants.h"
#include "../../headers/initfuncs.h"
#include "../../headers/misc.h"
#include "../../headers/geometry.h"
#include <math.h>

double rho_0(int cell[DIM_NUM-1]) {
	if( 1.5 * cell[0] - x1cellnum / 3 <= cell[1]) {
		return 100.;
	} else {
		return 1.;
	}
}


double p_0(int cell[DIM_NUM-1]) {
	if(1.5 * cell[0] - x1cellnum / 3 <= cell[1]) {
		return 133.3333333333;
	} else {
		return 1.;
	}
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
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
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
