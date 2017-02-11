#include <math.h>
#include "../headers/constants.h"
#include "../headers/misc.h"
#include "../headers/geometry.h"
#include "../headers/timestep.h"
#include "../headers/update.h"
#include <stdlib.h>
#include <stdio.h>

double t_now = t_init;
double h(double *Pstate_local, double x0, int cell[3]) {
	return 1. + gamma_const / (gamma_const - 1.) * Pstate_local[1] / Pstate_local[0];
}


double v_i(double *Pstate_local, int index, double x0, int cell[3]) {
	// contravariant components

	double *beta;
	double v_i_return;

	beta = shift(x0, cell);

	v_i_return = Pstate_local[2+index] / (lapse(x0, cell) * Pstate_local[2]) + 
		beta[index - 1] / lapse(x0, cell);
	
	free(beta);

	return v_i_return;
}


double W(double *Pstate_local, double x0, int cell[3]) {
	return lapse(x0, cell) * Pstate_local[2]; 
}


double D(double *Pstate_local, double x0, int cell[3]) {
	return Pstate_local[0] * W(Pstate_local, x0, cell);
}


double E(double *Pstate_local, double x0, int cell[3]) {
	return Pstate_local[0] * h(Pstate_local, x0, cell) * 
		pow(W(Pstate_local, x0, cell), 2.) - Pstate_local[1];
}


void S_k(double *Pstate_local, double result[3], double x0, int cell[3]) {
	double v[3];

	for(int i = 1; i < 4; i++) {
		v[i-1] = v_i(Pstate_local, i, x0, cell);
	}

	double h_here = h(Pstate_local, x0, cell);
	double W2 = pow(W(Pstate_local, x0, cell), 2.);

	for(int i = 1; i < 4; i++) {
		result[i-1] = Pstate_local[0] * h_here * W2 * lower_index_space(v, i, x0, cell);
	}
}


double T_mu_nu(double* Pstate_local, int dir1, int dir2, double x0, int cell[3]) {
	return Pstate_local[0] * h(Pstate_local, x0, cell) * Pstate_local[dir1 + 2] * Pstate_local[dir2 + 2] + 
		Pstate_local[1] * inverse_metric_full(dir1, dir2, x0, cell);
}


double* prim_to_cons_local(double *Pstate_local, double x0, int cell[3]) {
	double *Ustate_local;
	Ustate_local = malloc(sizeof(double) * 5);
	double S_i[3];
	
	S_k(Pstate_local, S_i, x0, cell);
	Ustate_local[0] = D(Pstate_local, x0, cell);
	Ustate_local[4] = E(Pstate_local, x0, cell) - Ustate_local[0];

	for(int i = 1; i < 4; i++) {
		Ustate_local[i] = S_i[i-1];
	}

	return Ustate_local;
}


void prim_to_cons(double *Ustate_ptr, double *Pstate_ptr, double x0) {

	int cell[3];
	double *Pstate_local;
	double S_i[3];
	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;
				
				Pstate_local = Pstate_ptr + P_offset(cell, 0);
				
				// define prim as P_ptr + P_offest(cell, 0)
				*(Ustate_ptr + U_offset(cell, 0)) = D(Pstate_local, x0, cell);
				*(Ustate_ptr + U_offset(cell, 4)) = E(Pstate_local, x0, cell) - *(Ustate_ptr + U_offset(cell, 0));

				S_k(Pstate_local, S_i, x0, cell);

				for(int index = 1; index < 4; index++) {
					*(Ustate_ptr + U_offset(cell, index)) = S_i[index-1];
				}
			}		
		}		
	}
}


double get_signals(double* Pstate_ptr, int plus, int dir, double x0, int cell[3]) {
	
	double *Pstate_left;
	int right_bias = 1;
	double* (*interpolation_method)(double*, int, int, double, int[3]);
	
	/*------------------------------------------------------------------------------------------------*/

	interpolation_method = piecewise_constant;

	/*------------------------------------------------------------------------------------------------*/

	Pstate_left = (*interpolation_method)(Pstate_ptr, right_bias, dir, x0, cell);

	double v[3];
	for(int i = 1; i < 4; i++) {
		v[i-1] = v_i(Pstate_left, i, x0, cell);
	}
	double vnormsq = inner_product_space(v, v, x0, cell);
	double p = Pstate_left[1];
	double rho = Pstate_left[0];
	double h_here = h(Pstate_left, x0, cell);
	double c = sqrt(gamma_const * p / (rho * h_here));
	double gii = inverse_metric_space(dir, dir, x0, cell);
	double a = lapse(x0, cell);

	cell[dir-1] += 1;

	double *Pstate_right;
	right_bias = 0;

	Pstate_right = (*interpolation_method)(Pstate_ptr, right_bias, dir, x0, cell);

	double vp1[3];
	for(int i = 1; i < 4; i++) {
		vp1[i-1] = v_i(Pstate_right, i, x0, cell);
	}
	double vnormsqp1 = inner_product_space(vp1, vp1, x0, cell);
	double pp1 = Pstate_right[1];
	double rhop1 = Pstate_right[0];
	double hp1 = h(Pstate_right, x0, cell);
	double cp1 = sqrt(gamma_const * pp1 / (rhop1 * hp1));
	double giip1 = inverse_metric_space(dir, dir, x0, cell);
	double ap1 = lapse(x0, cell);
 

	cell[dir-1] -= 1;



	double part = c * sqrt((1. - vnormsq) * (gii * (1. - vnormsq * pow(c, 2.)) - 
		pow(v[dir-1], 2.) * (1. - pow(c, 2.))));

	double partp1 = cp1 * sqrt((1. - vnormsqp1) * (giip1 * (1. - vnormsqp1 * pow(cp1, 2.)) - 
		pow(vp1[dir-1], 2.) * (1. - pow(cp1, 2.))));

	double *shiftvect;
	shiftvect = shift(x0, cell);

	double full, fullp1;
	
	if(plus) {
		full = a / (1. - vnormsq * pow(c, 2.0)) * (v[dir-1] * (1. - pow(c, 2.0)) + part) - shiftvect[dir-1];
		fullp1 = ap1 / (1. - vnormsqp1 * pow(cp1, 2.0)) * (vp1[dir-1] * (1. - pow(cp1, 2.0)) + partp1) - shiftvect[dir-1];
	} else {
		full = a / (1. - vnormsq * pow(c, 2.0)) * (v[dir-1] * (1. - pow(c, 2.0)) - part) - shiftvect[dir-1];
		fullp1 = ap1 / (1. - vnormsqp1 * pow(cp1, 2.0)) * (vp1[dir-1] * (1. - pow(cp1, 2.0)) - partp1) - shiftvect[dir-1];
	}
	free(Pstate_left);
	free(Pstate_right);
	free(shiftvect);

	if(fabs(full) > fabs(fullp1)) {
		return full;
	} else {
		return fullp1;
	}
}


double *dx_min = NULL;



double get_dt_max(double *Pstate_ptr, double x0, double *dxmin) {
	int cell[3];
	double c_max = 0.;
	double temp_plus;
	double temp_minus;
	if(dxmin == NULL) {
		dxmin = malloc(sizeof(double));
		*dxmin = calculate_dx_min(t_init);
	}

	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;

				for(int dir = 1; dir < 4; dir++) {
					temp_plus = fabs(get_signals(Pstate_ptr, 1, dir, x0, cell));
					temp_minus = fabs(get_signals(Pstate_ptr, 0, dir, x0, cell));

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

	// FIX THIS TRY HARMONIC MEAN

	return cfl_number * *dxmin / c_max;
}


double* build_F_i(double *Pstate_local, int dir, double x0, int cell[3]) {
	double *F_i_return;
	F_i_return = malloc(sizeof(double) * 5);

	double *beta;
	beta = shift(x0, cell);
	double S_i[3];
	S_k(Pstate_local, S_i, x0, cell);

	double alpha = lapse(x0, cell);
	double factor = v_i(Pstate_local, dir, x0, cell) - beta[dir-1] / alpha;
	F_i_return[0] = D(Pstate_local, x0, cell) * factor;
	F_i_return[1] = S_i[0] * factor;
	F_i_return[2] = S_i[1] * factor;
	F_i_return[3] = S_i[2] * factor;
	F_i_return[4] = (E(Pstate_local, x0, cell) - D(Pstate_local, x0, cell)) * 
		factor + Pstate_local[1] * v_i(Pstate_local, dir, x0, cell);

	F_i_return[dir] += Pstate_local[1];

	free(beta);

	return F_i_return;
}


double* piecewise_constant(double *Pstate_ptr, int right_bias, int dir, double x0, int cell[3]) {
	double *Pstate_return;
	Pstate_return = malloc(sizeof(double) * 6);
	for(int comp = 0; comp < 6; comp++) {
		Pstate_return[comp] = *(Pstate_ptr + P_offset(cell, comp));
	}

	return Pstate_return;
}


double* rhlle(double *Pstate_ptr, int dir, double x0, int cell[3]) {
	
	double *Pstate_left;
	double *Pstate_right;

	double* (*interpolation_method)(double*, int, int, double, int[3]);

	/*------------------------------------------------------------------------------------------------*/

	interpolation_method = piecewise_constant;

	/*------------------------------------------------------------------------------------------------*/

	Pstate_left = (*interpolation_method)(Pstate_ptr, 1, dir, x0, cell);
	
	cell[dir-1] += 1;

	Pstate_right = (*interpolation_method)(Pstate_ptr, 0, dir, x0, cell);
	
	cell[dir-1] -= 1;
	

	double *Uleft;
	double *Uright;
	double *Fleft;
	double *Fright;
	double *Freturn;
	double alpha_plus = get_signals(Pstate_ptr, 1, dir, x0, cell);
	double alpha_minus = get_signals(Pstate_ptr, 0, dir, x0, cell);
	
	Freturn = malloc(sizeof(double) * 5);
	Fleft = build_F_i(Pstate_left, dir, x0, cell);
	Fright = build_F_i(Pstate_right, dir, x0, cell);
	Uleft = prim_to_cons_local(Pstate_left, x0, cell);
	Uright = prim_to_cons_local(Pstate_right, x0, cell);

	if(alpha_plus <= 0.) {
		alpha_plus = 0.;
	}

	if(alpha_minus >= 0.) {
		alpha_minus = 0.;
	}

	for(int comp = 0; comp < 5; comp++) {
		Freturn[comp] = (alpha_plus * Fleft[comp] - alpha_minus * Fright[comp] - 
			alpha_plus * alpha_minus * (Uleft[comp] - Uright[comp])) / (alpha_plus - alpha_minus);
		Freturn[comp] *= sqrt(-metric_det(x0, cell));
	}

	free(Uleft);
	free(Uright);
	free(Fleft);
	free(Fright);
	free(Pstate_left);
	free(Pstate_right);

	return Freturn;
}


double* fluxes(double *Pstate_ptr, double x0, int cell[3]) {
	double* flux_sum;
	double* flux_left;
	double* flux_right;
	flux_sum = malloc(sizeof(double) * 5);
	
	for(int comp = 0; comp < 5; comp++) {
		flux_sum[comp] = 0.;
	}

	
	for(int dir = 1; dir < 4; dir++) {
	
		cell[dir-1] -= 1;
		flux_left = rhlle(Pstate_ptr, dir, x0, cell);
		cell[dir-1] += 1;
		flux_right = rhlle(Pstate_ptr, dir, x0, cell);
	
		for(int comp = 0; comp < 5; comp++) {
			flux_sum[comp] += (flux_left[comp] - flux_right[comp]) / dx_i(dir, x0, cell);
		}

		free(flux_left);
		free(flux_right);

	}

	return flux_sum;
}


double* sources(double *Pstate_ptr, double x0, int cell[3]) {
	double *source_return;
	source_return = malloc(sizeof(double) * 5);
	double christoffel_sum = 0.;
	double Tsum1 = 0.;
	double Tsum2 = 0.;
	double *Pstate_local;
	Pstate_local = Pstate_ptr + P_offset(cell, 0);

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
				source_return[comp] += T_mu_nu(Pstate_local, dir1, dir2, x0, cell) * 
					(metric_derivative(dir2, comp, dir1, x0, cell) - christoffel_sum);
				christoffel_sum = 0.;
			}
		}
	}
	for(int dir1 = 0; dir1 < 4; dir1++) {
		Tsum1 += T_mu_nu(Pstate_local, dir1, 0, x0, cell) * log_lapse_derivative(dir1, x0, cell);

		for(int dir2 = 0; dir2 < 4; dir2++) {
			Tsum2 += T_mu_nu(Pstate_local, dir1, dir2, x0, cell) * christoffel_symbol(0, dir1, dir2, x0, cell);
		}
	}

	source_return[4] = lapse(x0, cell) * (Tsum1 - Tsum2);

	for(int comp = 0; comp < 5; comp++) {
		source_return[comp] *= sqrt(-metric_det(x0, cell));
	}
	return source_return;
}


double update(double *Ustate_ptr, double *Pstate_ptr, double x0) {
	double dt = get_dt_max(Pstate_ptr, x0, dx_min);
	timestep_basic(Ustate_ptr, Pstate_ptr, dt, x0);
	return dt;
}


void pressureeq(double pbar, double *Ustate_ptr, double result[2], double x0, int cell[3]) {
	double *Ustate_local;
	Ustate_local = Ustate_ptr + U_offset(cell, 0);

	double D = Ustate_local[0];
	double tau = Ustate_local[4];
	double v_i[3];
	double v_upper_i[3];
	double v_sum = 0.;
	
	for(int i = 1; i < 4; i++) {
		v_i[i-1] = Ustate_local[i] / (tau + D + pbar);
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
	double e = (tau + D * (1. - W) + pbar * (1. - pow(W, 2.))) / (D * W);

	result[0] = pbar - (gamma_const - 1.) * rho * e;
	result[1] = v_sum * gamma_const * pbar / (rho + gamma_const * pbar / (gamma_const - 1.)) - 1.;

}


void bracket(double p, double result[3], double *Ustate_ptr, double x0, int cell[3]) {
	double upperp[2];
	double lowerp[2];
	double n = 2.0;


	while(1) {
	
		pressureeq(p * n, Ustate_ptr, upperp, x0, cell);
	
		pressureeq(p / n, Ustate_ptr, lowerp, x0, cell);

		if(sign(lowerp[0]) != sign(upperp[0])) {
		
			result[1] = p * n;
			result[0] = p / n;
		
			result[2] = sign(upperp[0]);
			break;

		} else {
			n *= 2.0;
		}

		if(n > 1.0e9) {
			printf("bracketing failed\n");
			break;
		}
	}

}


double psolve(double *Ustate_ptr, double *Pstate_ptr, double x0, int cell[3]) {
	double *Pstate_local;
	Pstate_local = Pstate_ptr + P_offset(cell, 0);
	int iterations = 0;
	double tolerance = 1e-8;
	double error = 1.0;
	double p = Pstate_local[1];
	double pfunc[2];
	double bounds[3];
	double delta;
	double *Ustate_local;
	
	bracket(p, bounds, Ustate_ptr, x0, cell);

	while(error >= tolerance) {
		pressureeq(p, Ustate_ptr, pfunc, x0, cell);
		delta = pfunc[0] / pfunc[1];
		if(fabs(delta) <= bounds[1] - bounds[0]) {
			error = fabs(delta);
		} else {
			error = bounds[1] - bounds[0];
		}

		if(p - delta < bounds[0] || p - delta > bounds[1] || isnan(delta)) {
			pressureeq((bounds[1] + bounds[0]) / 2., Ustate_ptr, pfunc, x0, cell);
			if(sign(pfunc[0]) != bounds[2]) {
				bounds[0] = (bounds[1] + bounds[0]) / 2.;
			} else {
				bounds[1] = (bounds[1] + bounds[0]) / 2.;
			}

			p = (bounds[1] + bounds[0]) / 2.;
		} else {
			p -= delta;
		}

		iterations += 1;
		if(iterations > 1000) {
			printf("psolve failed");
			break;
		}
	}

	return p;
}


void cons_to_prim(double *Ustate_ptr, double *Pstate_ptr, double x0) {
	int cell[3];
	double *Ustate_local;
	double D, tau, p, W, v_sum, u0;
	double v_i[3];
	double v_i_upper[3];
	double *shift_vect;
	
	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;

				v_sum = 0.;
			
				Ustate_local = Ustate_ptr + U_offset(cell, 0);

				p = psolve(Ustate_ptr, Pstate_ptr, x0, cell);
				D = Ustate_local[0];
				tau = Ustate_local[4];

				for(int i = 1; i < 4; i++) {
					v_i[i-1] = Ustate_local[i] / (tau + D + p);
				}
				
				for(int i = 1; i < 4; i++) {
					v_i_upper[i-1] = 0.;
					for(int j = 1; j < 4; j++){
						v_i_upper[i-1] += inverse_metric_full(i, j, x0, cell) * v_i[j-1];
					}
					v_sum += v_i[i-1] * v_i_upper[i-1];
				}

				W = pow(1. - v_sum, -.5);
				u0 = W / lapse(x0, cell);
				shift_vect = shift(x0, cell);

				*(Pstate_ptr + P_offset(cell, 0)) = D / W;
				*(Pstate_ptr + P_offset(cell, 1)) = p;
				*(Pstate_ptr + P_offset(cell, 2)) = u0;
				*(Pstate_ptr + P_offset(cell, 3)) = u0 * (lapse(x0, cell) * v_i_upper[0] - shift_vect[0]);
				*(Pstate_ptr + P_offset(cell, 4)) = u0 * (lapse(x0, cell) * v_i_upper[1] - shift_vect[1]);
				*(Pstate_ptr + P_offset(cell, 5)) = u0 * (lapse(x0, cell) * v_i_upper[2] - shift_vect[2]);
				
				free(shift_vect);
			}		
		}		
	}
}

