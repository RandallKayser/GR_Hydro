#include <math.h>
#include "update.h"
#include "geometry.h"
#include "timestep.h"

double h(double *Pstate_local, double x0, int cell[DIM_NUM-1]) {
	return 1. + gamma / (gamma - 1.) * *(Pstate_local + 1) / *Pstate_local;
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


double T_mu_nu(double* Pstate_local, int dir1, int dir2, double x0, int cell[DIM_NUM-1]) {
	return *Pstate_local * h(Pstate_local, x0, cell) * *(Pstate_local + dir1 + 2) * (Pstate_local + dir2 + 2) + 
		*(Pstate_local + 1) * inverse_metric_full(dir1, dir2, x0, cell);
}


double* prim_to_cons_local(double *Pstate_local, double x0, int cell[DIM_NUM-1]) {
	double *Ustate_local;
	Ustate_local = malloc(sizeof(double) * 5);
	Ustate_local[0] = D(Pstate_local, x0, cell);
	Ustate_local[1] = E(Pstate_local, x0, cell) - Ustate_local[0];
	for(int i = 1; i < DIM_NUM, i++) {
		Ustate_local[1 + i] = S_k(Pstate_local, i, x0, cell);
	}

	return Ustate_local;
}


void prim_to_cons(double *Ustate_ptr, double *Pstate_ptr, double x0) {

	int cell[DIM_NUM-1];
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
				// define prim as P_ptr + P_offest(cell, 0)
				*(Ustate_ptr + U_offset(cell, 0)) = D(Pstate_local, x0, cell);
				*(Ustate_ptr + U_offset(cell, 1)) = E(Pstate_local, x0, cell) - *(Ustate_ptr + U_offset(cell, 0));

				for(int index = 1; index < DIM_NUM; index++) {
					*(Ustate_ptr + U_offset(cell, 1 + index)) = S_k(Pstate_local, index, x0, cell);
				}
			}		
		}		
	}
}


double get_signals(double* Pstate_ptr, int plus, int dir, double x0, int cell[DIM_NUM-1]) {
	
	double Pstate_left[5];
	int right_bias = 1;
	double* (*interpolation_method)(double*, int, int, double, int[DIM_NUM-1]);
	
/*------------------------------------------------------------------------------------------------*/

	interpolation_method = piecewise_constant;

/*------------------------------------------------------------------------------------------------*/

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


double* rhlle(double *Pstate_ptr, int dir, double x0, int cell[DIM_NUM-1]) {
	
	double Pstate_left[5];
	double Pstate_right[5];
	double* (*interpolation_method)(double*, int, int, double, int[DIM_NUM-1]);

/*------------------------------------------------------------------------------------------------*/

	interpolation_method = piecewise_constant;

/*------------------------------------------------------------------------------------------------*/

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


double* fluxes(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	double* flux_sum;
	double* flux_left;
	double* flux_right;
	flux_sum = malloc(sizeof(double) * 5);


	for(int comp = 0; comp < 5; comp++) {
		flux_sum[comp] = 0.;
		
		for(int dir = 1; dir < 4; dir++) {
			cell[dir-1] -= 1;
			flux_left = rhlle(Pstate_ptr, dir, x0, cell);
			cell[dir-1] += 1;
			flux_right = rhlle(Pstate_ptr, dir, x0, cell);

			flux_sum[comp] += (flux_left[comp] - flux_right[comp]) / dx_i(dir, cell);
			// look at this use of dx_i -- might cause issues

			free(flux_left);
			free(flux_right);
	
		}
	}

	return flux_sum;
}


double* sources(double *Pstate_ptr, double x0, int cell[DIM_NUM-1]) {
	double *source_return;
	source_return = malloc(sizeof(double) * 5);
	double christoffel_sum = 0.;
	double Tsum1 = 0.;
	double Tsum2 = 0.;
	for(int comp = 0; comp < 5; comp++) {
		source_return[comp] = 0.;
	}
	for(int comp = 1; comp < 4; comp++) {
		for(int dir1 = 0; dir1 < 4; dir1++) {
			for(int dir2 = 0; dir2 < 4; dir2++) {
				for(int delta = 0; delta < 4; delta++) {
					christoffel_sum += christoffel_symbol(delta, dir1, dir2, x0, cell) * 
						metric_full(delta, comp, x0, cell);
				}
				source_return[comp] += T_mu_nu(Pstate_ptr, dir1, dir2, x0, cell) * 
					(metric_derivative(dir2, comp, dir1) - christoffel_sum);
				christoffel_sum = 0.;
			}
		}
	}
	for(int dir1 = 0; dir1 < 4; dir1++) {
		Tsum1 += T_mu_nu(dir1, 0, x0, cell) * log_lapse_derivative(dir1, x0, cell);

		for(int dir2 = 0; dir2 < 4; dir2++) {
			Tsum2 += T_mu_nu(dir1, dir2, x0, cell) * christoffel_symbol(0, dir1, dir2, x0, cell);
		}
	}

	source_return[4] = lapse(x0, cell) * (Tsum1 - Tsum2);
}



double* piecewise_constant(*Pstate_ptr, int right_bias, int dir, double x0, int cell[DIM_NUM-1]) {
	double *Pstate_return;
	Pstate_return = malloc(sizeof(double) * 5);
	for(comp = 0; comp < 5; comp++) {
		Pstate_return[comp] = *(Pstate_ptr + P_offset(cell, comp));
	}

	return Pstate_return;
}


void update(double *Ustate_ptr, double *Pstate_ptr, double x0) {
	timestep_basic(Ustate_ptr, Pstate_ptr, get_dt_max(Pstate_ptr, x0));

}


void pressureeq(double pbar, double *Ustate_ptr, double result[2], double x0, int cell[DIM_NUM-1]) {
	double D = *Ustate_ptr;
	double tau = *(Ustate_ptr + 1);
	double v_i[3];
	double v_upper_i[3];
	double v_sum = 0.;
	for(int i = 1; i < 4; i++) {
		v_i[i-1] = *(Ustate_ptr + i) / (tau + D + pbar);
	}
	for(int i = 1; i < 4; i++) {
		v_upper_i[i-1] = 0.;
		for(int j = 1; j < 4; j++){
			v_upper_i[i-1] += inverse_metric_full(i, j, x0, cell) * v_i[j-1];
		}
		v_sum += v_i[i-1] * v_upper_i[i-1];
	}

	double W = pow(1. - v_sum, -.5);
	double rho = D / W;
	double e = (tau + pbar + D * (1. - W) + pbar * (1. - pow(W, 2.))) / (D * W);

	result[0] = pbar - (gamma - 1.) * rho * e;
	result[1] = v_sum * gamma * pbar / (rho + gamma * pbar / (gamma - 1.));

}


double bisection(double pguess, double *Ustate_ptr, double x0, int cell[DIM_NUM-1]) {
	int flag = 1;
	double n = .1;
	double fu[2];
	double fm[2];
	double fl[2];
	double upperbound;
	double lowerbound;
	double midpoint;
	double tolerance = 1e-5;

	while(flag) {
		pressureeq(pguess + n, Ustate_ptr, fu, x0, cell);
		pressureeq(pguess - n, Ustate_ptr, fl, x0, cell);
		if(sign(fu[0]) != sign(fl[0])) {
			upperbound = n;
			lowerbound = 0;
			flag = 0;

		} else {
			n += 1.0;
			if(n >= 1000.0) {
				printf("%f, %f, %f\n", D, S, tau);
				flag = 0;
			}
			free(fu);
			free(fl);
		}
	}

	while(upperbound - lowerbound >= tolerance) {

		midpoint = (upperbound + lowerbound) / 2.0;
		fm = pressureeq(midpoint, D, S, tau);

		if(sign(fu[0]) != sign(fm[0])) {
			lowerbound = midpoint;
			free(fl);
			free(fm);
			fl = pressureeq(lowerbound, D, S, tau);
		} else {
			upperbound = midpoint;
			free(fu);
			free(fm);
			fu = pressureeq(upperbound, D, S, tau);
		}
	}

	return midpoint;
}


void psolve(double Ustate[100][3][xcellnum], int timestep, int position, double primitives[100][xcellnum][4]) {
	double D = Ustate[timestep][0][position];
	double S = Ustate[timestep][1][position];
	double tau = Ustate[timestep][2][position];
	int iterations = 0;
	double tolerance = 1e-10;
	double error = 1.0;
	double p = primitives[(timestep+99) % 100][position][2];
	double* pfunc;

	p = bisection(p, D, S, tau, position);

	while(error >= tolerance && iterations < 100) {
		pfunc = pressureeq(p, D, S, tau);
		error = fabs(pfunc[0] / pfunc[1]);
		p -= pfunc[0] / pfunc[1];
		iterations += 1;
	}

	double v = S / (tau + D + p);
	double rho = D * sqrt(1.0 - v * v);
	double h = 1.0 + agamma * p / ((agamma - 1.0) * rho);

	primitives[timestep][position][0] = rho;
	primitives[timestep][position][1] = v;
	primitives[timestep][position][2] = p;
	primitives[timestep][position][3] = h;

}
























