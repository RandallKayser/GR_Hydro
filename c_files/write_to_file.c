#include "../headers/constants.h"
#include "../headers/misc.h"
#include "../headers/write_to_file.h"
#include <stdlib.h>

void P_dump(FILE *f, double *Pstate_ptr, double x0) {
	int cell[3];
	double *Pstate_local;

	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;
				
				Pstate_local = Pstate_ptr + P_offset(cell, 0);
				fprintf(f, "%i, %i, %i, %f, %f, %f, %f, %f, %f, %f\n", cell[0], cell[1], cell[2], x0,
					Pstate_local[0], Pstate_local[1], Pstate_local[2], Pstate_local[3], 
					Pstate_local[4], Pstate_local[5]);

			}	
		}	
	}
}