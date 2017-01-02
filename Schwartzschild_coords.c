#include "constants.h"
#include "geometry.h"
#include <math.h>

double mass = 1;
double eps = 1e-6;
double r_min = 2. * mass + eps;
double r_max;

// Standard Swartzschild coordinates, {t, r, theta, phi}
// logarithmic cells in r
// delta_theta, delta_phi constant 

double schwarzschild(int dir1, int dir2, double position[4]) {
	if(dir1 != dir2) {
		return 0;
	} else if(dir1 == 0) {
		return -1. + 2. * mass / position[1];
	} else if(dir1 == 1) {
		return pow(1. - 2. * mass / position[1], -1.);
	} else if(dir3 == 2) {
		return pow(position[1], 2.);
	} else if(dir1 == 3) {
		return pow(position[1] * sin(position[2]), 2.);
	}
}

double get_position(int dir, int cell[DIM_NUM-1]) {
	if(dir == 1) {
		return r_min * pow(r_max / r_min, cell[0] / x1cellnum);
	} else if(dir == 2) {
		return 2 * M_PI * cell[1] / x2cellnum;
	} else if(dir == 3) {
		return M_PI * cell[2] / x3cellnum;
	} 
}

double dx_i(int dir, int cell[DIM_NUM-1]) {
	if(dir == 1) {
		double xp1;
		double xm1;

		cell[0] += 1;
		xp1 = get_position(1, cell);
		cell[0] -= 2;
		xm1 = get_position(1, cell);
		cell[0] += 1;
		return (xp1 - xm1) / 2.;
	} else if(dir == 2) {
		return 2 * M_PI / x2cellnum * get_position(1, cell) * sin(get_position(3, cell));
	} else if(dir == 3) {
		return M_PI / x3cellnum * get_position(1, cell);
	}
}

double calculate_dx_min() {
	dx_min = 1e30;
	double cell[DIM_NUM-1];
	double dx_temp;
	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;

				for(dir = 1; dir < DIM_NUM; dir++) {
					dx_temp = dx_i(dir, cell);

					if(dx_temp < dx_min) {
						dx_min = dx_temp
					}

				} 
			}
		}
	}

	return dx_min;
}


const static double dx_min = calculate_dx_min();


double metric_full(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	double position[4];

	position[0] = x0;
	position[1] = get_position(1, cell);
	position[2] = get_position(2, cell);
	position[3] = get_position(3, cell);

	return schwarzschild(dir1, dir2, position);
}


double inner_product(double a[DIM_NUM], double b[DIM_NUM], double x0, int cell[DIM_NUM-1]) {

	double sum = 0.;

	for(int dir = 0; dir < DIM_NUM; dir++) {
		sum += metric_full(dir, dir, x0, cell) * a[dir] * b[dir];
	}

	return sum;
}


double metric_det(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	double r = get_position(1, cell);
	double theta = get_position(2, cell);

	return pow(pow(r, 2.) * sin(theta), 2.);
}


double metric_space(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	return metric_full(int dir1, int dir2, double x0, int cell[DIM_NUM-1]);
}


double inverse_metric_space(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	if(dir1 != dir2) {
		return 0.;
	} else {
		return pow(schwarzschild(dir1, dir2, x0, cell), -1.);
	}
}


double lapse(double x0, int cell[DIM_NUM-1]) {
	return pow(-metric_full(0, 0, x0, cell[DIM_NUM-1]), .5);
}


double *shift(double x0, int cell[DIM_NUM-1]) {
	return 0.;
}


double inner_product_space(double a[DIM_NUM-1], double b[DIM_NUM-1], double x0, int cell[DIM_NUM-1]) {
	double sum = 0.;

	for(int dir = 1; dir < DIM_NUM; dir++) {
		sum += metric_space(dir, dir, x0, cell) * a[dir-1] * b[dir-1];
	}

	return sum;
}	


double lower_index_space(double v[DIM_NUM-1], int index) {
	double sum = 0.;

	for(int dir = 1; dir < DIM_NUM; dir++) {
		sum += metric_space(index, dir, x0, cell) * v[dir-1];
	}

	return sum;
}	

double volume(double x0, int cell[DIM_NUM-1]) {
	double x1 = get_position(1, x0, cell);
	double x3 = get_position(3, x0, cell);
	double dx1 = dx_i(1, cell);
	double dx2 = dx_i(2, cell);
	double dx3 = dx_i(3, cell);

	double volume = (pow(x1 + dx1, 3.) + pow(x1, 3.)) / 3. * dx2 * (cos(x3) - cos(x3 + dx3));

	return volume;
}




