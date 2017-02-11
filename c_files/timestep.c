#include "../headers/timestep.h"
#include "../headers/update.h"
#include "../headers/boundary.h"
#include "../headers/misc.h"
#include <stdlib.h>

void timestep_basic(double *Ustate_ptr, double *Pstate_ptr, double dt, double x0) {
	
	int cell[DIM_NUM-1];
	double *flux_vect;
	double *source_vect;

	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;
				flux_vect = fluxes(Pstate_ptr, x0, cell);
				source_vect = sources(Pstate_ptr, x0, cell);

				for(int comp = 0; comp < DIM_NUM + 1; comp++) {
					*(Ustate_ptr + U_offset(cell, comp)) += dt * (flux_vect[comp]); //+ source_vect[comp]);
				}

				enforce_bc(Ustate_ptr, Pstate_ptr, x0, x1type, x2type, x3type);

				free(flux_vect);
				free(source_vect);
			}
		}
	}
}