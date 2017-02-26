#include "../headers/constants.h"
#include "../headers/initfuncs.h"
#include "../headers/misc.h"
#include "../headers/timestep.h"
#include "../headers/update.h"
#include "../headers/write_to_file.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void print_diagnose(double *Ustate_ptr, double *Pstate_ptr, double x0) {
	int cell[3];
	double *Ustate_local;
	double *Pstate_local;	
	for(int i = 0; i < 2 * ghost_num + x1cellnum; i++) {
		for(int j = 0; j < 2 * ghost_num + x2cellnum; j++) {
			for(int k = 0; k < 2 * ghost_num + x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;

				Ustate_local = Ustate_ptr + U_offset(cell, 0);
				Pstate_local = Pstate_ptr + P_offset(cell, 0);

				printf("%f   %i %i %i, %f %f %f %f %f, %f %f %f %f %f %f\n", x0, cell[0], cell[1], cell[2], 
					Ustate_local[0], Ustate_local[1], Ustate_local[2], Ustate_local[3], Ustate_local[4],
					Pstate_local[0], Pstate_local[1], Pstate_local[2], Pstate_local[3], Pstate_local[4], Pstate_local[5]);
			}	
		}		
	}
}

int main(void) {
	double *U_ptr;
	double *U1_ptr;
	double *P_ptr;
	double t_now = t_init;
	double *t_now_ptr = &t_now;
	double dump_interval = (t_final - t_init) / 500.;
	double checkpoint = 0.;

	if(ghost_num == 1) {
		interpolation_method = piecewise_constant;
		RK_type = forward_euler;
		r_solver = rhlle;
	} else if(ghost_num == 2) {
		r_solver = rhllc;
		interpolation_method = plm;
		RK_type = RK2_TVD;
	} else if(ghost_num == 3) {
		interpolation_method = ppm;
		RK_type = RK3_TVD;
	}

	FILE *f = fopen("output/first_test.txt", "w");

	if(f == NULL) {
		printf("NULL ERROR");
	}

	U_ptr = malloc(sizeof(double) * 5 * (2 * ghost_num + x1cellnum) * (2 * ghost_num + x2cellnum) * (2 * ghost_num + x3cellnum));
	U1_ptr = malloc(sizeof(double) * 5 * (2 * ghost_num + x1cellnum) * (2 * ghost_num + x2cellnum) * (2 * ghost_num + x3cellnum));

	P_ptr = malloc(sizeof(double) * 6 * (2 * ghost_num + x1cellnum) * (2 * ghost_num + x2cellnum) * (2 * ghost_num + x3cellnum));

	Pinit(P_ptr);
	while(*t_now_ptr < t_final) {
//		if(*t_now_ptr >= t_init + checkpoint * dump_interval) {
			P_dump(f, P_ptr, *t_now_ptr);
//			checkpoint += 1.;
//		}
//		print_diagnose(U_ptr, P_ptr, *t_now_ptr);
		(*RK_type)(U_ptr, U1_ptr, P_ptr, t_now_ptr);
		printf("time elapsed = %f\n", *t_now_ptr);
	}

	return 0;
}
