#include "constants.h"
#include "initfuncs.h"
#include "misc.h"
#include "geometry.h"
#include "timestep.h"
#include "update.h"
#include <stdlib.h>
#include <stdio.h>

int main(void) {
	double t_now = t_init;
	U_ptr = malloc(sizeof(double) * 5 * x0cellnum * x1cellnum * x2cellnum * x3cellnum);
	P_ptr = malloc(sizeof(double) * 6 * x0cellnum * x1cellnum * x2cellnum * x3cellnum);

	Pinit(P_ptr);
	printf("dxmin = %f\n", calculate_dx_min(t_init));
	while(t_now < t_final) {
		update(U_ptr, P_ptr, t_now);
		printf("%f\n", t_now);
	}
/* above, keep
------------------------------------------------------------------------------------------------------------------
   below, test funcs */
	int cell[3];
	
	cell[0] = 14;
	cell[1] = 5;
	cell[2] = 7;

	double *rhlle_test;
	double *flux_test;
	double *source_test;

	rhlle_test = rhlle(P_ptr, 1, t_init, cell);
	flux_test = fluxes(P_ptr, t_init, cell);
	source_test = sources(P_ptr, t_init, cell);

	printf("%f %f %f %f \n", get_position(1, t_init, cell), get_position(2, t_init, cell), get_position(3, t_init, cell));
	printf("%f %f %f \n", dx_i(1, t_init, cell), dx_i(2, t_init, cell), dx_i(3, t_init, cell));
	printf("%f %f %f %f %f\n", rhlle_test[0], rhlle_test[1], rhlle_test[2], rhlle_test[3], rhlle_test[4]);
	printf("%f %f %f %f %f\n", flux_test[0], flux_test[1], flux_test[2], flux_test[3], flux_test[4]);
	printf("%f %f %f %f %f\n", source_test[0], source_test[1], source_test[2], source_test[3], source_test[4]);

	return 0;
}