#include <math.h>
#include "update.h"
#include "geometry.h"

double h(double *Pstate_local, double x0, int cell[DIM_NUM-1]) {
	return 1. + (gamma - 1.) / gamma * *(Pstate_local + 1) / *Pstate_local;
}


double v_i(double *Pstate_local, int index, double x0, int cell[DIM_NUM-1]) {
	

	double *beta;
	double v_i_return;

	beta = shift(x0, cell);
	v_i_return = *(Pstate_local + 2 + index) / (lapse(x0, cell) * 
		*(Pstate_local + 2) + beta[index - 1] / lapse(x0, cell);
	
	free(beta);

	return v_i_return;
}


double W(double *Pstate_local, double x0, int cell[DIM_NUM-1]) {
	return lapse(x0, cell) * *(Pstate_local + 2); 
}


double D(double *Pstate_local, double x0, int cell[DIM_NUM-1]) {

	return *Pstate_local * W(Pstate_local, x0, cell);
}


double E(double *Pstate_local, double x0, int cell[DIM_NUM-1]) {
	
	return *Pstate_local * h(Pstate_local, x0, cell) * 
		pow(W(Pstate_local, x0, cell), 2.) - *(Pstate_local + 1);
}


double S_k(double *Pstate_local, int index, double x0, int cell[DIM_NUM-1]) {
	double v[DIM_NUM-1];

	for(int i = 1; i < DIM_NUM, i++) {
		v[i-1] = v_i(Pstate_local, i, x0, cell);
	}

	return *Pstate_local * h(Pstate_local, x0, cell) * pow(W(Pstate_local, x0, cell), 2.) * 
		lower_index(v, index)
}


/* old version, probably to be deleted

double h(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	return 1. + (gamma - 1.) / gamma * *(Pstate_ptr + P_offset(cell, 1)) / 
		*(Pstate_ptr + P_offset(cell, 0));
}


double v_i(double *Pstate_ptr, int index, double x0, int cell[DIM_NUM-1]) {
	

	double *beta;
	double v_i_return;

	beta = shift(x0, cell);
	v_i_return = *(Pstate_ptr + P_offset(cell, 2+index)) / (lapse(x0, cell) * 
		*(Pstate_ptr + P_offset(cell, 2)) + beta[index - 1] / lapse(x0, cell);
	
	free(beta);

	return v_i_return;
}


double W(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	
	double v[DIM_NUM-1];
	double W_return;

	for(int dir = 1; dir < DIM_NUM; dir++) {
		v[dir-1] = v_i(Pstate_ptr, dir, x0, cell);
	}

	return pow(1. - inner_3_product(v, v), -.5); 
}


double D(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {

	return *(Pstate_ptr + P_offset(cell, 0)) * W(Pstate_ptr, x0, cell);
}


double E(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	
	return *(Pstate_ptr + P_offset(cell, 0)) * h(Pstate, x0, cell) * 
		pow(W(Pstate, x0, cell), 2.) - *(Pstate_ptr + P_offset(cell, 1));
}


double S_k(double *Pstate_ptr, int index, double x0, int cell[DIM_NUM-1]) {
	double v[DIM_NUM-1];

	for(int i = 1; i < DIM_NUM, i++) {
		v[i-1] = v_i(Pstate_ptr, i, x0, cell);
	}

	return *(Pstate_ptr + P_offset(cell, 0)) * h(Pstate_ptr, x0, cell) * pow(W(Pstate_ptr, x0, cell), 2.) * 
		lower_index(v, index)
}
*/


double* prim_to_cons_local(double *Pstate_local, double x0, int cell[DIM_NUM-1]) {
	double *Ustate_local;
	Ustate_local = malloc(sizeof(double) * 5);
	Ustate_local[0] = D(Pstate_local, x0, cell);
	Ustate_local[1] = E(Pstate_local, x0, cell);
	for(int i = 1; i < DIM_NUM, i++) {
		Ustate_local[1 + i] = S_k(Pstate_local, i, x0, cell);
	}

	return Ustate_local;
}


void prim_to_cons(double *Ustate_ptr, double *Pstate_ptr, double x0) {

	double cell[DIM_NUM-1];
	double Pstate_local[DIM_NUM+2];

	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;
				for(int comp = 0; comp < DIM_NUM+2; comp++) {
					Pstate_local[comp] = *(Pstate_ptr + P_offset(cell, comp));
				}

				*(Ustate_ptr + U_offset(cell, 0)) = D(Pstate_local, x0, cell);
				*(Ustate_ptr + U_offset(cell, 1)) = E(Pstate_local, x0, cell);

				for(int index = 1; index < DIM_NUM; index++) {
					*(Ustate_ptr + U_offset(cell, 1 + index)) = S_k(Pstate_local, index, x0, cell);
				}
			}		
		}		
	}
}


double get_signals(double* Pstate_ptr, int plus, int dir, double x0, 
	int cell[DIM_NUM-1], double* (*interpolation_method)(double*, int, int, double, int[DIM_NUM-1]))) {
	
	double Pstate_left[5];
	int right_bias = 1;

	Pstate_left = (*interpolation_method)(Pstate_ptr, right_bias, dir, x0, cell);

	double v[DIM_NUM-1];
	for(int i = 1; i < DIM_NUM-1; i++) {
		v[i-1] = v_i(Pstate_left, i, x0, cell)
	}
	double vnormsq = inner_product_space(v, v, x0, cell);
	double p = *(Pstate_left + 1);
	double rho = *Pstate_left;
	double h = h(Pstate_left, x0, cell);
	double c = sqrt(gamma * p / (rho * h));
	double gii = inverse_metric_space(dir, dir, x0, cell);
	double a = lapse(x0, cell);

	cell[dir-1] += 1;

	double Pstate_right[5];
	right_bias = 0;

	Pstate_right = (*interpolation_method)(Pstate_ptr, right_bias, dir, x0, cell);

	double vp1[DIM_NUM-1];
	for(int i = 1; i < DIM_NUM-1; i++) {
		vp1[i-1] = v_i(Pstate_ptr, i, x0, cell)
	}
	double vnormpsq1 = inner_product_space(vp1, vp1, x0, cell);
	double pp1 = *(Pstate_right + 1);
	double rhop1 = *Pstate_right;
	double hp1 = h(Pstate_right, x0, cell);
	double cp1 = sqrt(gamma * pp1 / (rhop1 * hp1));
	double giip1 = inverse_metric_space(dir, dir, x0, cell);
	double ap1 = lapse(x0, cell);

	cell[dir-1] -= 1;

	double vbar = (v[dir-1] + vp1[dir-1]) / 2.;
	double cbar = (c + cp1) / 2.;
	double vnormsqbar = (vnormsq + vnormsqp1) / 2.;
	double giibar = (gii + giip1) / 2.
	double abar = (a + ap1) / 2.;

	double part = cbar * sqrt((1. - vnormsqbar) * (giibar * (1. - vnormsqbar * pow(cbar, 2.)) - 
		pow(vbar, 2.) * (1. - pow(cbar, 2.))));
	double *shiftvect;
	shiftvect = shift(x0, cell);

	if(plus) {
		double full = abar / (1. - vnormsqbar * pow(cbar, 2.0)) * (vbar * (1. - pow(cbar, 2.0)) + part) - shiftvect[dir-1];
	} else {
		double full = abar / (1. - vnormsqbar * pow(cbar, 2.0)) * (vbar * (1. - pow(cbar, 2.0)) - part) - shiftvect[dir-1];
	}

	free(shiftvect);

	return full;
}


double get_dt_max(double *Pstate_ptr, double x0) {
	double cell[DIM_NUM-1];
	double c_max = 0;
	double temp_plus;
	double temp_minus;
	double temp_dx;
	for(int i = 0; i < x1cellnum; i++) {
		for(int j = 0; j < x2cellnum; j++) {
			for(int k = 0; k < x3cellnum; k++) {
				cell[0] = i;
				cell[1] = j;
				cell[2] = k;

				for(dir = 1; dir < DIM_NUM; dir++) {
					temp_plus = fabs(get_signals(Pstate_ptr, 1, dir, x0, cell));
					temp_minus = fabs(get_signals(Pstate_ptr, 0, dir, x0, cell));
					temp_dx = dx_i(dir, cell);

					if(temp_plus > c_max) {
						c_max = temp_plus;
					}

					if(temp_minus > c_max) {
						c_max = temp_minus;
					}

				} 
			}
		}
	}

	return cfl_number * dx_min / c_max;
}


double* build_F_i(double *Pstate_local, int dir, double x0, int cell[DIM_NUM-1]) {
	double *F_i_return;
	F_i_return = malloc(sizeof(double) * 5);

	double *beta;
	beta = shift(x0, cell);
	
	double alpha = lapse(x0, cell);
	double factor = v_i(Pstate_local, dir, x0, cell) - beta[dir-1] / alpha;
	F_i_return[0] = D(Pstate_local, x0, cell) * factor;
	F_i_return[1] = S_k(Pstate_local, 1, x0, cell) * factor;
	F_i_return[2] = S_k(Pstate_local, 2, x0, cell) * factor;
	F_i_return[3] = S_k(Pstate_local, 3, x0, cell) * factor;
	F_i_return[4] = (E(Pstate_local, x0, cell) - D(Pstate_local, x0, cell)) * 
		factor + *(Pstate_local + 1) * v_i(Pstate_local, dir, x0, cell);

	F_i_return[dir] += *(Pstate_local + 1);

	free(beta);

	return F_i_return;
}


double* rhlle(double *Pstate_ptr, int dir, double* 
	(*interpolation_method)(double*, int, int, double, int[DIM_NUM-1]), double x0, int cell[DIM_NUM-1]) {
	
	double Pstate_left[5];
	double Pstate_right[5];
	Pstate_left = (*interpolation_method)(Pstate_ptr, 1, dir, x0, cell);
	cell[dir-1] += 1;

	Pstate_right = (*interpolation_method)(Pstate_ptr, 0, dir, x0, cell);
	cell[dir-1] -= 1;
	
	// likely going to limit convergence -- metric not evaluated at edges inside get_signals()
	double Uleft[5];
	double Uright[5];
	double Fleft[5];
	double Fright[5];
	double Freturn[5];
	double alpha_plus = get_signals(Pstate_ptr, 1, dir, x0, cell, interpolation_method);
	double alpha_minus = get_signals(Pstate_ptr, 0, dir, x0, cell, interpolation_method);
	Fleft = build_F_i(Pstate_left, dir, x0, cell);
	Fright = build_F_i(Pstate_right, dir, x0, cell);
	Uleft = prims_to_cons_local(Pstate_left, x0, cell);
	Uright = prims_to_cons_local(Pstate_right, x0, cell);

	if(alpha_plus <= 0.) {
		alpha_plus = 0.;
	}

	if(alpha_minus >= 0.) {
		alpha_minus = 0.;
	}

	for(int comp = 0; comp < 5; comp++) {
		Freturn[comp] = (alpha_plus * Fleft[comp] - alpha_minus * Fright[comp] - 
			alpha_plus * alpha_minus * (Uleft[comp] - Uright[comp]) / (alpha_plus - alpha_minus);
	}

	free(Fleft);
	free(Fright);
	free(Pstate_left);
	free(Pstate_right);

	return Freturn;
}





































