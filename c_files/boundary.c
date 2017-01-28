#include "../headers/constants.h"
#include "../headers/boundary.h"

void enforce_bc(double *Ustate_ptr, int x1type, int x2type, int x3type) {
	int cell[3];
	void (*bc)();

	if(x1type == 0) {
		bc = bc_periodic;
	} else if(x1type == 1) {
		bc = bc_fixed_value;
	} else if(x1type == 2) {
		bc = bc_zero_gradient;
	} else if(x1type == 3) {
		bc = bc_reflecting;
	}
	for(int i = 0; i < ghost_num; i++) {
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				
			}
		}
	}
}