#include "timestep.h"
#include "update.h"
#include "misc.h"

void timestep_basic(double *Ustate_ptr, double *Pstate_ptr, double dt, double x0) {
	
	int cell[DIM_NUM-1];
	double *flux_vect;
	double *source_vect;

	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;
				flux_vect = fluxes(Pstate_ptr, x0, cell);
				source_vect = sources(Pstate_ptr, x0, cell);

				for(int comp = 0; comp < DIM_NUM + 1; comp++) {
					*(Ustate_ptr + U_offset(cell, comp)) += dt * (flux_vect[comp] + source_vect[comp]);
				}
			}
		}
	}
}


