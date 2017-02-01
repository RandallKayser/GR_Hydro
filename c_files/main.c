#include "../headers/constants.h"
#include "../headers/initfuncs.h"
#include "../headers/misc.h"
#include "../headers/geometry.h"
#include "../headers/timestep.h"
#include "../headers/update.h"
#include "../headers/boundary.h"
#include "../headers/write_to_file.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int main(void) {
	
	double t_now = t_init;
	FILE *f;
	printf("WE GOT HERE1");
	f = fopen("output/first_test.txt", "w");
	printf("WE GOT HERE");
	U_ptr = malloc(sizeof(double) * 5 * (2 * ghost_num + x1cellnum) * (2 * ghost_num + x2cellnum) * (2 * ghost_num + x3cellnum));
	P_ptr = malloc(sizeof(double) * 6 * (2 * ghost_num + x1cellnum) * (2 * ghost_num + x2cellnum) * (2 * ghost_num + x3cellnum));

	Pinit(P_ptr);

	while(t_now < t_final) {
		P_dump(f, P_ptr, t_now);
		prim_to_cons(U_ptr, P_ptr, t_now);
		enforce_bc(U_ptr, P_ptr, t_now, x1type, x2type, x3type);
		t_now += update(U_ptr, P_ptr, t_now);
		cons_to_prim(U_ptr, P_ptr, t_now);
		printf("time elapsed = %f\n", t_now);
	}

	return 0;
}

