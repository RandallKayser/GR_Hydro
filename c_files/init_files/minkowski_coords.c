#include "../../headers/geometry.h"
#include <stdio.h>
#include <stdlib.h>

const static double length_x = 10.;
const static double length_y = 10.;
const static double length_z = 10.;
const static double dx = length_x / ((double) x1cellnum - 1.);
const static double dy = length_y / ((double) x2cellnum - 1.);
const static double dz = length_z / ((double) x3cellnum - 1.);

double get_position(int dir, double x0, int cell[DIM_NUM-1]) {
	if(dir == 1) {
		return (double) cell[0] * dx + dx/2.;
	} else if(dir == 2) {
		return (double) cell[1] * dy + dy/2.;
	} else if(dir == 3) {
		return (double) cell[2] * dz + dz/2.;
	} else {
		printf("index error in get_position");
		return 0.;
	}
}


double dx_i(int dir, double x0, int cell[DIM_NUM-1]) {
	if(dir == 1) {
		return dx;
	} if(dir == 2) {
		return dy;
	} if(dir == 3) {
		return dz;
	} else {
		printf("index error in dx_i");
		return 0.;
	}
}


double calculate_dx_min(double x0) {
	if(dx <= dy) {
		if(dx <= dz) {
			return dx;
		} else {
			return dz;
		}
	} else if(dy <= dz) {
		return dy;
	} else {
		return dz;
	}
}


double metric_full(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	if(dir1 != dir2) {
		return 0.;
	}
	if(dir1 == 0) {
		return -1.;
	} else {
		return 1.;
	}
}


double metric_derivative(int dir1, int dir2, int dir3, double x0, int cell[DIM_NUM-1]){
	// dir1 dir2 are the metric indices, dir3 the partial derivative index
	return 0.;
}


double christoffel_symbol(int dir1, int dir2, int dir3, double x0, int cell[DIM_NUM-1]) {
	// dir1 as top index
	return 0.;	
}


double inverse_metric_full(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	return metric_full(dir1, dir2, x0, cell);
}


double inner_product(double a[DIM_NUM], double b[DIM_NUM], double x0, int cell[DIM_NUM-1]) {
	double sum = 0.;
	for(int dir = 0; dir < DIM_NUM; dir++) {
		sum += a[dir] * b[dir] * metric_full(dir, dir, x0, cell);
	}
	return sum;
}


double metric_det(double x0, int cell[DIM_NUM-1]) {
	return -1.;
}

double metric_space(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	if(dir1 == 0 || dir2 == 0) {
		printf("index error in metric_space, looking for time coordinate.\n");
		return 0.;
	} else {
		return metric_full(dir1, dir2, x0 ,cell);
	}
}


double inverse_metric_space(int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	return metric_space(dir1, dir2, x0, cell);
}


double lapse(double x0, int cell[DIM_NUM-1]) {
	return 1.;
}


double log_lapse_derivative(int dir, double x0, int cell[DIM_NUM-1]) {
	return 0.;
}


double *shift(double x0, int cell[DIM_NUM-1]) {
	double *beta_return;
	beta_return = malloc(sizeof(double) * 3);
	beta_return[0] = 0.;
	beta_return[1] = 0.;
	beta_return[2] = 0.;

	return beta_return;
}


double inner_product_space(double a[DIM_NUM-1], double b[DIM_NUM-1], double x0, int cell[DIM_NUM-1]) {
	double sum = 0.;
	for(int dir = 1; dir < DIM_NUM; dir++) {
		sum += a[dir-1] * b[dir-1] * metric_full(dir, dir, x0, cell);
	}
	return sum;
}

double lower_index_space(double v[DIM_NUM-1], int index, double x0, int cell[DIM_NUM-1]) {
	return v[index-1];
}

double volume(double x0, int cell[DIM_NUM-1]) {
	return dx * dy * dz;
}



