#include <math.h>
#include "update.h"
#include "geometry.h"

double h(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	return 1. + (gamma - 1.) / gamma * *(Pstate_ptr + P_offset(cell, 1)) / 
		*(Pstate_ptr + P_offset(cell, 0));
}


double v_i(double *Pstate_ptr, int index, double x0, int cell[DIM_NUM-1]) {
	

	double *beta;
	double v_i_return;

	beta = shift(x0, cell);
	v_i_return = *(Pstate_ptr + P_offset(cell, 2+index)) / (lapse(x0, cell) * 
		*(Pstate_ptr + P_offset(cell, 2)) + beta[index - 1] / lapse(x0, cell);
	
	free(beta);

	return v_i_return;
}


double W(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	
	double v[3];
	double W_return;

	for(int dir = 1; dir < 4; dir++) {
		v[dir-1] = v_i(Pstate_ptr, dir, x0, cell);
	}

	return pow(1. - inner_3_product(v, v), -.5); 
}


double D(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {

	return *(Pstate_ptr + P_offset(cell, 0)) * W(Pstate_ptr, x0, cell);
}


double E(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	
	return *(Pstate_ptr + P_offset(cell, 0)) * h(Pstate, x0, cell) * 
		pow(W(Pstate, x0, cell), 2.) - *(Pstate_ptr + P_offset(cell, 1));
}


double S_k(double *Pstate_ptr, int index, double x0, int cell[DIM_NUM-1]) {
	double v[DIM_NUM-1];

	for(int i = 1; i < DIM_NUM, i++) {
		v[i-1] = v_i(Pstate_ptr, i, x0, cell);
	}

	return *(Pstate_ptr + P_offset(cell, 0)) * h(Pstate_ptr, x0, cell) * pow(W(Pstate_ptr, x0, cell), 2.) * 
		lower_index(v, index)
}


void prim_to_cons(double *Ustate_ptr, double *Pstate_ptr, double x0) {

	double cell[DIM_NUM-1];

	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;

				*(Ustate_ptr + U_offset(cell, 0)) = D(Pstate_ptr, x0, cell);
				*(Ustate_ptr + U_offset(cell, 1)) = E(Pstate_ptr, x0, cell);

				for(int index = 1; index < DIM_NUM; index++) {
					*(Ustate_ptr + U_offset(cell, 1 + index)) = S_k(Pstate_ptr, index, x0, cell);
				}
			}		
		}		
	}
}


double get_signals(int plus, double *Pstate_ptr, int dir, double x0, int cell[DIM_NUM-1]) {
	if(plus) {
		if(dir == 1) {
			return blah1plus;
		} else if(dir == 2) {
			return blah2plus;
		} else if(dir == 3) {
			return blah3plus;
		}
	} else {
		if(dir == 1) {
			return blah1minus;
		} else if(dir == 2) {
			return blah2minus;
		} else if(dir == 3) {
			return blah3minus;
		}
	}
	
}










































