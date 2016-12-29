#include "timestep.h"

void timestep_basic(double *Ustate_ptr, double *Pstate_ptr, double *Ustate_temp, 
	double dt, void (*subroutine)(double*, double*, double* double, double[DIM_NUM-1])) {
	
	int cell[DIM_NUM-1];
	
	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;
				
				for(int comp = 0; comp < DIM_NUM + 1; comp++) {
					*(Ustate_temp + U_offset(cell, comp)) = *(Ustate_ptr + U_offset(cell, comp)) + 
						dt * (*subroutine)(Ustate_ptr, Pstate_ptr, Ustate_temp, dt, cell);
					/* subroutines: add_fluxes, add_sources, */
				}
			}
		}
	}
}


