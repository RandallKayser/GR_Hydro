#include "headers/constants.h"
#include "headers/geometry.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Standard Swartzschild coordinates, {t, r, theta, phi}
// logarithmic cells in r
// delta_theta, delta_phi constant 

double schwarzschild(int dir1, int dir2, double position[4]) {
	if(dir1 != dir2) {
		return 0.;
	} else if(dir1 == 0) {
		return -1. + 2. * mass / position[1];
	} else if(dir1 == 1) {
		return pow(1. - 2. * mass / position[1], -1.);
	} else if(dir1 == 2) {
		return pow(position[1], 2.);
	} else if(dir1 == 3) {
		return pow(position[1] * sin(position[2]), 2.);
	} else {
		printf("\n\n\nmetric index error\n\n\n");
		return 0.;
	}
}


double get_position(int dir, double x0, int cell[DIM_NUM-1]) {
	if(dir == 1) {
		return r_min * pow(r_max / r_min, (double) cell[0] / x1cellnum);
	} else if(dir == 2) {
		return M_PI * cell[1] / x2cellnum;
	} else if(dir == 3) {
		return 2 * M_PI * cell[2] / x3cellnum;
	} else {
		printf("\n\n\nget_position_index_error\n\n\n");
		return 0.;
	}
}


double dx_i(int dir, double x0, int cell[DIM_NUM-1]) {
	if(dir == 1) {
		double xp1;
		double xm1;

		cell[0] += 1;
		xp1 = get_position(1, x0, cell);
		cell[0] -= 2;
		xm1 = get_position(1, x0, cell);
		cell[0] += 1;
		return (xp1 - xm1) / 2.;
	} else if(dir == 2) {
		return 2 * M_PI / x2cellnum * get_position(1, x0, cell) * sin(get_position(2, x0, cell));
	} else if(dir == 3) {
		return M_PI / x3cellnum * get_position(1, x0, cell);
	} else {
		printf("\n\n\ndx_i index error\n\n\n");
		return 0.;
	}
}


double calculate_dx_min(double x0) {
	double dx_min = 1e30;
	int cell[DIM_NUM-1];
	double dx_temp;
	for(int i = 1; i < x1cellnum - 1; i++) {
		for(int j = 1; j < x2cellnum - 1; j++) {
			for(int k = 1; k < x3cellnum - 1; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;



				for(int dir = 1; dir < DIM_NUM; dir++) {
					if(dx_i(dir, x0, cell) == 0.) {
						printf("ERROR AT %i %i %i\n", cell[0], cell[1], cell[2]);
					}
					dx_temp = dx_i(dir, x0, cell);

					if(dx_temp < dx_min) {
						dx_min = dx_temp;
					}

				} 
			}
		}
	}

	return dx_min;
}


double metric_full(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	double position[4];

	position[0] = x0;
	position[1] = get_position(1, x0, cell);
	position[2] = get_position(2, x0, cell);
	position[3] = get_position(3, x0, cell);

	return schwarzschild(dir1, dir2, position);
}


double metric_derivative(int dir1, int dir2, int dir3, double x0, int cell[DIM_NUM-1]) {
	double r = get_position(1, x0, cell);
	double theta = get_position(2, x0, cell);

	if(dir3 == 0 || dir3 == 3) {
		return 0.;
	} else if(dir1 != dir2) {
		return 0.;
	} else if(dir1 == 0 && dir3 == 1) {
		return - 2. * mass / pow(r, 2.);
	} else if(dir1 == 0 && dir3 == 2) {
		return 0.;
	} else if(dir1 == 1 && dir3 == 1) {
		return 2 * mass / pow(r - 2 * mass, 2.);
	} else if(dir1 == 1 && dir3 == 2) {
		return 0.;
	} else if(dir1 == 2 && dir3 == 1) {
		return 2 * r;
	} else if(dir1 == 2 && dir3 == 2) {
		return 0;
	} else if(dir1 == 3 && dir3 == 1) {
		return 2 * r * pow(sin(theta), 2.);
	} else if(dir1 == 3 && dir3 == 2) {
		return pow(r, 2.) * 2 * sin(theta) * cos(theta);
	} else {
		printf("\n\n Error in metric_derivative\n\n");
		return 0.;
	}
}


double christoffel_symbol(int dir1, int dir2, int dir3, double x0, int cell[DIM_NUM-1]) {
	// dir1 as top index
	double r = get_position(1, x0, cell);
	double theta = get_position(2, x0, cell);

	if(dir2 < dir3) {
		int temp = dir3;
		dir3 = dir2;
		dir2 = temp;
	}

	if(dir1 == 0) {
		if(dir2 == 1 && dir3 == 0) {
			return mass / (pow(r, 2.) - 2. * mass * r);
		}
	} else if(dir1 == 1) {
		if(dir2 == dir3) {
			if(dir2 == 0) {
				return mass / pow(r, 2.) * (1. - 2. * mass / r);
			} else if(dir2 == 1) {
				return mass / (pow(r, 2.) * (1. - 2. * mass / r));
			} else if(dir2 == 2) {
				return -r * (1. - 2. * mass / r);
			} else if(dir2 == 3) {
				return -r * pow(sin(theta), 2.) * (1. - 2. * mass / r);
			}
		}
	} else if(dir1 == 2) {
		if(dir2 == 2 && dir3 == 1) {
			return pow(r, -1.);
		} else if(dir2 == 3 && dir3 == 3) {
			return -sin(theta) * cos(theta);
		}
	} else if(dir1 == 3) {
		if(dir2 == 3 && dir3 == 1) {
			return pow(r, -1.);
		} else if(dir2 == 3 && dir3 == 2) {
			return cos(theta) / sin(theta);
		}
	}

	return 0.;
}


double inverse_metric_full(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	if(dir1 != dir2) {
		return 0.;
	} else {
		return pow(metric_full(dir1, dir2, x0, cell), -1.);
	}
}


double inner_product(double a[DIM_NUM], double b[DIM_NUM], double x0, int cell[DIM_NUM-1]) {
	// takes 2 contravariant vectors
	double sum = 0.;

	for(int dir = 0; dir < DIM_NUM; dir++) {
		sum += metric_full(dir, dir, x0, cell) * a[dir] * b[dir];
	}

	return sum;
}


double metric_det(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	double r = get_position(1, x0, cell);
	double theta = get_position(2, x0, cell);

	return pow(pow(r, 2.) * sin(theta), 2.);
}


double metric_space(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	return metric_full(dir1, dir2, x0, cell);
}


double inverse_metric_space(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	if(dir1 != dir2) {
		return 0.;
	} else {
		return pow(metric_space(dir1, dir2, x0, cell), -1.);
	}
}


double lapse(double x0, int cell[DIM_NUM-1]) {
	return pow(-metric_full(0, 0, x0, cell), .5);
}


double log_lapse_derivative(int dir, double x0, int cell[DIM_NUM-1]) {
	if(dir == 1) {
		return mass * pow(lapse(x0, cell) * get_position(1, x0, cell), 2.);
	} else {
		return 0.;
	}
}


double *shift(double x0, int cell[DIM_NUM-1]) {
	double *shift_vect;
	shift_vect = malloc(sizeof(double) * 4);
	for(int i = 0; i < DIM_NUM; i++) {
		shift_vect[i] = 0.;
	}

	return shift_vect;
}


double inner_product_space(double a[DIM_NUM-1], double b[DIM_NUM-1], double x0, int cell[DIM_NUM-1]) {
	double sum = 0.;

	for(int dir = 1; dir < DIM_NUM; dir++) {
		sum += metric_space(dir, dir, x0, cell) * a[dir-1] * b[dir-1];
	}

	return sum;
}	


double lower_index_space(double v[DIM_NUM-1], int index, double x0, int cell[DIM_NUM-1]) {
	double sum = 0.;

	for(int dir = 1; dir < DIM_NUM; dir++) {
		sum += metric_space(index, dir, x0, cell) * v[dir-1];
	}

	return sum;
}	

double volume(double x0, int cell[DIM_NUM-1]) {
	double x1 = get_position(1, x0, cell);
	double x3 = get_position(3, x0, cell);
	double dx1 = dx_i(1, x0, cell);
	double dx2 = dx_i(2, x0, cell);
	double dx3 = dx_i(3, x0, cell);

	double volume = (pow(x1 + dx1, 3.) + pow(x1, 3.)) / 3. * dx2 * (cos(x3) - cos(x3 + dx3));

	return volume;
}