#include "../headers/constants.h"
#include "../headers/boundary.h"
#include "../headers/update.h"
#include "../headers/misc.h"
#include <stdlib.h>
#include <stdio.h>

void bc_periodic(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]) {
	double *Pstate_edge;
	double *Ustate_edge;
	int plus = 0;

	if(cell[dir-1] < ghost_num) {
		cell[dir-1] += dircellnum;
		plus = 1;
	} else {
		cell[dir-1] -= dircellnum;
	}
	
	Pstate_edge = Pstate_ptr + P_offset(cell, 0);

	Ustate_edge = prim_to_cons_local(Pstate_edge, x0, cell);

	if(plus == 0) {
		cell[dir-1] += dircellnum;
	} else {
		cell[dir-1] -= dircellnum;
	}

	for(int comp = 0; comp < 1 + DIM_NUM; comp++) {
		*(Ustate_ptr + U_offset(cell, comp)) = Ustate_edge[comp];
		*(Pstate_ptr + P_offset(cell, comp)) = Pstate_edge[comp];
	}
	*(Pstate_ptr + P_offset(cell, 5)) = Pstate_edge[5];

	free(Ustate_edge);

}


void bc_fixed_value(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]) {

}

void bc_zero_gradient(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]) {
	double *Pstate_edge;
	double *Ustate_edge;

	int temp = cell[dir-1];
	
	if(cell[dir-1] < dircellnum) {
		cell[dir-1] = ghost_num;
	} else {
		cell[dir-1] = dircellnum + ghost_num - 1;
	}

	Pstate_edge = Pstate_ptr + P_offset(cell, 0);

	Ustate_edge = prim_to_cons_local(Pstate_edge, x0, cell);

	cell[dir-1] = temp;

	for(int comp = 0; comp < 5; comp++) {
		*(Ustate_ptr + U_offset(cell, comp)) = Ustate_edge[comp];
		*(Pstate_ptr + P_offset(cell, comp)) = Pstate_edge[comp];
	}

	*(Pstate_ptr + P_offset(cell, 5)) = Pstate_edge[5];

	free(Ustate_edge);
}


void bc_reflecting(double *Ustate_ptr, double *Pstate_ptr, int dir, int dircellnum, double x0, int cell[3]) {

}

void enforce_bc(double *Ustate_ptr, double *Pstate_ptr, double x0, int x1type, int x2type, int x3type) {
	int cell[3];
	void (*bc)(double*, double*, int, int, double, int[3]);

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
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;

				(*bc)(Ustate_ptr, Pstate_ptr, 1, x1cellnum, x0, cell);
			}
		}

		cell[0] = x1cellnum + 2 * ghost_num - (i + 1);
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;

				(*bc)(Ustate_ptr, Pstate_ptr, 1, x1cellnum, x0, cell);
			}
		}
	}

	if(x2type == 0) {
		bc = bc_periodic;
	} else if(x2type == 1) {
		bc = bc_fixed_value;
	} else if(x2type == 2) {
		bc = bc_zero_gradient;
	} else if(x2type == 3) {
		bc = bc_reflecting;
	}

	for(int i = 0; i < ghost_num; i++) {	
		cell[1] = i;
		for(int j = ghost_num; j < x1cellnum + ghost_num; j++) {
			cell[0] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;

				(*bc)(Ustate_ptr, Pstate_ptr, 2, x2cellnum, x0, cell);
			}
		}

		cell[1] = x2cellnum + 2 * ghost_num - (i + 1);
		for(int j = ghost_num; j < x1cellnum + ghost_num; j++) {
			cell[0] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;

				(*bc)(Ustate_ptr, Pstate_ptr, 2, x2cellnum, x0, cell);
			}
		}
	}

	if(x3type == 0) {
		bc = bc_periodic;
	} else if(x3type == 1) {
		bc = bc_fixed_value;
	} else if(x3type == 2) {
		bc = bc_zero_gradient;
	} else if(x3type == 3) {
		bc = bc_reflecting;
	}

	for(int i = 0; i < ghost_num; i++) {
		cell[2] = i;
		for(int j = ghost_num; j < x1cellnum + ghost_num; j++) {
			cell[0] = j;
			for(int k = ghost_num; k < x2cellnum + ghost_num; k++) {
				cell[1] = k;

				(*bc)(Ustate_ptr, Pstate_ptr, 3, x3cellnum, x0, cell);		
			}
		}

		cell[2] = x3cellnum + 2*ghost_num - (i + 1);		
		for(int j = ghost_num; j < x1cellnum + ghost_num; j++) {
			cell[0] = j;
			for(int k = ghost_num; k < x2cellnum + ghost_num; k++) {
				cell[1] = k;

				(*bc)(Ustate_ptr, Pstate_ptr, 3, x3cellnum, x0, cell);		
			}
		}
	}
}

