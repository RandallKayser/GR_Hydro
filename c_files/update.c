#include <math.h>
#include "../headers/constants.h"
#include "../headers/misc.h"
#include "../headers/geometry.h"
#include "../headers/update.h"
#include <stdlib.h>
#include <stdio.h>


double h(double *Pstate_local, double x0, int cell[3]) {
	return 1.0 + gamma_const * Pstate_local[1] / ((gamma_const - 1.0) * Pstate_local[0]);
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
		pow(W(Pstate_local, x0, cell), 2.0) - Pstate_local[1];
}


void S_k(double *Pstate_local, double result[3], double x0, int cell[3]) {
	double v[3];

	for(int i = 1; i < 4; i++) {
		v[i-1] = v_i(Pstate_local, i, x0, cell);
	}

	double h_here = h(Pstate_local, x0, cell);
	double W2 = pow(W(Pstate_local, x0, cell), 2.0);

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

	Pstate_left = (*interpolation_method)(Pstate_ptr, 0, dir, x0, cell);	

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


	double *Pstate_right;

	Pstate_right = (*interpolation_method)(Pstate_ptr, 1, dir, x0, cell);
	
	cell[dir-1] += 1;
	
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


	double part = c * sqrt((1.0 - vnormsq) * (gii * (1.0 - vnormsq * pow(c, 2.0)) - 
		pow(v[dir-1], 2.0) * (1.0 - pow(c, 2.0))));

	double partp1 = cp1 * sqrt((1.0 - vnormsqp1) * (giip1 * (1.0 - vnormsqp1 * pow(cp1, 2.0)) - 
		pow(vp1[dir-1], 2.0) * (1.0 - pow(cp1, 2.0))));

	cell[dir-1] -= 1;
	
	double *shiftvect;
	shiftvect = shift(x0, cell);

	double full, fullp1;
	
	if(plus) {
		full = a / (1.0 - vnormsq * pow(c, 2.0)) * (v[dir-1] * (1.0 - pow(c, 2.0)) + part) - shiftvect[dir-1];
		fullp1 = ap1 / (1.0 - vnormsqp1 * pow(cp1, 2.0)) * (vp1[dir-1] * (1.0 - pow(cp1, 2.0)) + partp1) - shiftvect[dir-1];
	} else {
		full = a / (1.0 - vnormsq * pow(c, 2.0)) * (v[dir-1] * (1.0 - pow(c, 2.0)) - part) - shiftvect[dir-1];
		fullp1 = ap1 / (1.0 - vnormsqp1 * pow(cp1, 2.0)) * (vp1[dir-1] * (1.0 - pow(cp1, 2.0)) - partp1) - shiftvect[dir-1];
	}
	free(Pstate_left);
	free(Pstate_right);
	free(shiftvect);

	if(plus) {
		if(full > fullp1) {
			return full;
		} else {
			return fullp1;
		}	
	} else {
		if(full > fullp1) {
			return fullp1;
		} else {
			return full;
		}
	}
}


double get_dt_max(double *Pstate_ptr, double x0, double dxmin) {
	int cell[3];
	double c_max = 0.;
	double temp_plus;
	double temp_minus;

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

	return cfl_number * dxmin / c_max;
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

	for(int comp = 0; comp < 5; comp++) {
		F_i_return[comp] *= sqrt(-metric_det(x0, cell));
	}
	free(beta);

	return F_i_return;
}


double* piecewise_constant(double *Pstate_ptr, int right_bias, int dir, double x0, int cell[3]) {
	double *Pstate_return;
	Pstate_return = malloc(sizeof(double) * 6);
	
	if(right_bias) {
		cell[dir-1] += 1;
	}

	for(int comp = 0; comp < 6; comp++) {
		Pstate_return[comp] = *(Pstate_ptr + P_offset(cell, comp));
	}
	if(right_bias) {
		cell[dir-1] -= 1;
	}

	return Pstate_return;
}


double minmod_flux_limiter(double r) {
	if(r <= 0.0) {
		return 0.;
	}

	double a = flux_limiter_theta * r;
	double b = .5 + r / 2.;
	double c = flux_limiter_theta;

	if(c <= a) {
		if(c <= b) {
			return c;
		} else {
			return b;
		}
	} else {
		if(a <= b) {
			return a;
		} else {
			return b;
		}
	}
}


double* plm(double *Pstate_ptr, int right_bias, int dir, double x0, int cell[3]) {
	double *Pstate_return;
	Pstate_return = malloc(sizeof(double) * 6);
	double r[6];	
	double *Pm1;
	double *P;
	double *Pp1;
	double *Pp2;

	P = Pstate_ptr + P_offset(cell, 0);
	cell[dir-1] += 1;
	Pp1 = Pstate_ptr + P_offset(cell, 0);

	if(right_bias) {
		cell[dir-1] += 1;
		Pp2 = Pstate_ptr + P_offset(cell, 0);

		for(int comp = 0; comp < 6; comp++) {
			if(Pp2[comp] != Pp1[comp]) {
				r[comp] = (Pp1[comp] - P[comp]) / (Pp2[comp] - Pp1[comp]);					
			} else {
				r[comp] = 0.;
			}
			Pstate_return[comp] = Pp1[comp] - .5 * minmod_flux_limiter(r[comp]) * (Pp1[comp] - P[comp]);	
		}
		
		cell[dir-1] -= 2;
	
	} else {
		cell[dir-1] -= 2;
		Pm1 = Pstate_ptr + P_offset(cell, 0);

		for(int comp = 0; comp < 6; comp++) {
			if(Pp1[comp] != P[comp]) {
				r[comp] = (P[comp] - Pm1[comp]) / (Pp1[comp] - P[comp]);
			} else {
				r[comp] = 0.;
			}
			Pstate_return[comp] = P[comp] + .5 * minmod_flux_limiter(r[comp]) * (Pp1[comp] - P[comp]);
		}	

		cell[dir-1] += 1;

	}

	return Pstate_return;

}


double* ppm(double *Pstate_ptr, int right_bias, int dir, double x0, int cell[3]) {
	return NULL;
}

double* rhlle(double *Pstate_ptr, int dir, double x0, int cell[3]) {
	double *Pstate_left;
	double *Pstate_right;

	Pstate_left = (*interpolation_method)(Pstate_ptr, 0, dir, x0, cell);
	Pstate_right = (*interpolation_method)(Pstate_ptr, 1, dir, x0, cell);

	double *Uleft;
	double *Uright;
	double *Fleft;
	double *Fright;
	double *Freturn;
	double alpha_plus = get_signals(Pstate_ptr, 1, dir, x0, cell);
	double alpha_minus = get_signals(Pstate_ptr, 0, dir, x0, cell);
	
	Freturn = malloc(sizeof(double) * 5);
	Fleft = build_F_i(Pstate_left, dir, x0, cell);
	Uleft = prim_to_cons_local(Pstate_left, x0, cell);
	cell[dir-1] += 1;
	Fright = build_F_i(Pstate_right, dir, x0, cell);
	Uright = prim_to_cons_local(Pstate_right, x0, cell);
	cell[dir-1] -= 1;
	if(alpha_plus <= 0.0) {
		alpha_plus = 0.;
	}

	if(alpha_minus >= 0.0) {
		alpha_minus = 0.;
	}

	for(int comp = 0; comp < 5; comp++) {
		Freturn[comp] = (alpha_plus * Fleft[comp] - alpha_minus * Fright[comp] - 
			alpha_plus * alpha_minus * (Uleft[comp] - Uright[comp])) / (alpha_plus - alpha_minus);
	}

	free(Uleft);
	free(Uright);
	free(Fleft);
	free(Fright);
	free(Pstate_left);
	free(Pstate_right);

	return Freturn;

}


double *rhllc(double *Pstate_ptr, int dir, double x0, int cell[DIM_NUM-1]) {
	double *Pstate_left;
	double *Pstate_right;

	Pstate_left = (*interpolation_method)(Pstate_ptr, 0, dir, x0, cell);
	Pstate_right = (*interpolation_method)(Pstate_ptr, 1, dir, x0, cell);

	double *Uleft;
	double *Uright;
	double Uhlle[5];
	double *Fleft;
	double *Fright;
	double Fhlle[5];
	double *Freturn;
	double alpha_plus = get_signals(Pstate_ptr, 1, dir, x0, cell);
	double alpha_minus = get_signals(Pstate_ptr, 0, dir, x0, cell);
	double alpha_star = 0.;

	Fleft = build_F_i(Pstate_left, dir, x0, cell);
	Uleft = prim_to_cons_local(Pstate_left, x0, cell);
	cell[dir-1] += 1;
	Fright = build_F_i(Pstate_right, dir, x0, cell);
	Uright = prim_to_cons_local(Pstate_right, x0, cell);
	cell[dir-1] -= 1;

	if(alpha_plus <= 0) {
		Freturn = malloc(sizeof(double) * 5);

		for(int comp = 0; comp < 5; comp++) {
			Freturn[comp] = Fright[comp];
		}
	} else if(alpha_minus >= 0) {
		Freturn = malloc(sizeof(double) * 5);

		for(int comp = 0; comp < 5; comp++) {
			Freturn[comp] = Fleft[comp];
		}
	} else {
		Freturn = malloc(sizeof(double) * 5);

		for(int comp = 0; comp < 5; comp++) {
			Fhlle[comp] = (alpha_plus * Fleft[comp] - alpha_minus * Fright[comp] - alpha_plus * alpha_minus * (Uleft[comp] - Uright[comp])) / (alpha_plus - alpha_minus);
			Uhlle[comp] = (alpha_plus * Uright[comp] - alpha_minus * Uleft[comp] + Fleft[comp] - Fright[comp]) / (alpha_plus - alpha_minus);
		}

		if(Fhlle[4] + Fhlle[0] != 0.0) {
			alpha_star = (Uhlle[4] + Uhlle[0] + Fhlle[dir] - 
				sqrt(pow(Uhlle[4] + Uhlle[0] + Fhlle[dir], 2.0) - 4.0 * (Fhlle[4] + Fhlle[0]) * Uhlle[dir])) / (2.0 * (Fhlle[4] + Fhlle[0]));
		} 

		double p_star;
		double D_star;
		double S_star[3];
		double E_star;

		if(alpha_star >= 0) {
			// LEFT STATE CHOSEN
			double A = alpha_minus * (Uleft[4] + Uleft[0]) - Uleft[dir];
			double B = Uleft[dir] * (alpha_minus - Pstate_left[dir + 2] / Pstate_left[2]) - Pstate_left[1];

			p_star = (A * alpha_star - B) / (1.0 + alpha_minus * alpha_star);

			D_star = Uleft[0] * (alpha_minus - Pstate_left[2 + dir] / Pstate_left[2]) / (alpha_minus - alpha_star);

			for(int index = 1; index < 4; index++) {
				S_star[index-1] = Uleft[index] * (alpha_minus - Pstate_left[2+dir]/Pstate_left[2]) / (alpha_minus - alpha_star);
			}

			E_star = ((Uleft[4] + Uleft[0]) * (alpha_minus - Pstate_left[2+dir] / Pstate_left[2]) + p_star * alpha_star - Pstate_left[1] * Pstate_left[2+dir] / Pstate_left[2]) / (alpha_minus - alpha_star);

			Freturn[0] = D_star * alpha_star;
			Freturn[1] = S_star[0] * alpha_star;
			Freturn[2] = S_star[1] * alpha_star;
			Freturn[3] = S_star[2] * alpha_star;
			Freturn[4] = (E_star) * alpha_star + p_star * alpha_star;

			Freturn[dir] = (E_star + p_star) * pow(alpha_star, 2.0) + p_star;

		} else {
			// RIGHT STATE CHOSEN
			cell[dir-1] += 1;
			double A = alpha_plus * (Uright[4] + Uright[0]) - Uright[dir];
			double B = Uright[dir] * (alpha_plus - Pstate_right[dir + 2] / Pstate_right[2]) - Pstate_right[1];

			p_star = (A * alpha_star - B) / (1.0 + alpha_plus * alpha_star);
			
			D_star = Uright[0] * (alpha_plus - Pstate_right[2 + dir] / Pstate_right[2]) / (alpha_plus - alpha_star);

			for(int index = 1; index < 4; index++) {
				S_star[index-1] = Uright[index] * (alpha_plus - Pstate_right[dir+2] / Pstate_right[2]) / (alpha_plus - alpha_star);
			}

			E_star = ((Uright[4] + Uright[0]) * (alpha_plus - Pstate_right[2+dir] / Pstate_right[2]) + p_star * alpha_star - Pstate_right[1] * Pstate_right[2+dir] / Pstate_right[2]) / (alpha_plus - alpha_star);
			
			Freturn[0] = D_star * alpha_star;
			Freturn[1] = S_star[0] * alpha_star;
			Freturn[2] = S_star[1] * alpha_star;
			Freturn[3] = S_star[2] * alpha_star;
			Freturn[4] = (E_star) * alpha_star + p_star * alpha_star;

			Freturn[dir] = (E_star + p_star) * pow(alpha_star, 2.0) + p_star;

			cell[dir-1] -= 1;

		}

	}

	free(Uleft);
	free(Uright);
	free(Fleft);
	free(Fright);

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
		flux_left = (*r_solver)(Pstate_ptr, dir, x0, cell);
		cell[dir-1] += 1;
		flux_right = (*r_solver)(Pstate_ptr, dir, x0, cell);
	
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


void update(double *Ustate_ptr, double *U1state_ptr, double *Ustate_target, double *Pstate_ptr, double rkfac1, double rkfac2, double dt, double x0) {
	
	int cell[DIM_NUM-1];
	double *flux_vect;
	double *source_vect;

	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;
				flux_vect = fluxes(Pstate_ptr, x0, cell);
				source_vect = sources(Pstate_ptr, x0, cell);

				for(int comp = 0; comp < DIM_NUM + 1; comp++) {
					*(Ustate_target + U_offset(cell, comp)) = rkfac1 * *(Ustate_ptr + U_offset(cell, comp)) + 
					rkfac2 * (*(U1state_ptr + U_offset(cell, comp)) + dt * (flux_vect[comp] + source_vect[comp]));
				}

				free(flux_vect);
				free(source_vect);
			}
		}
	}
}


void pressureeq(double pbar, double *Ustate_local, double result[2], double x0, int cell[3]) {

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
	double W = pow(1.0 - v_sum, -.5);
	double rho = D / W;
	double e = (tau + D * (1.0 - W) + pbar * (1.0 - pow(W, 2.0))) / (D * W);

	result[0] = pbar - (gamma_const - 1.0) * rho * e;
	result[1] = v_sum * gamma_const * pbar / (rho + gamma_const * pbar / (gamma_const - 1.0)) - 1.;

}


void bracket(double p, double result[3], double *Ustate_local, double x0, int cell[3]) {
	double upperp[2];
	double lowerp[2];
	double n = 2.0;


	while(n <= 1.0e9) {
	
		pressureeq(p * n, Ustate_local, upperp, x0, cell);
	
		pressureeq(p / n, Ustate_local, lowerp, x0, cell);

		if(sign(lowerp[0]) != sign(upperp[0])) {
		
			result[1] = p * n;
			result[0] = p / n;
		
			result[2] = sign(upperp[0]);
			break;

		} else {
			n *= 2.0;
		}

	}

	if(n > 1.0e9) {
		printf("bracket failed at p = %.15f, D = %f, S = %f, tau = %.15f, x0 = %f, cell[0] = %i\n", p, Ustate_local[0], Ustate_local[1], Ustate_local[2], x0, cell[0]);
	}

}


double psolve(double *Ustate_local, double *Pstate_local, double x0, int cell[3]) {
	int iterations = 0;
	double tolerance = 1e-10;
	double error = 1.0;
	double p = Pstate_local[1];
	double pfunc[2];
	double bounds[3];
	double delta;

	bracket(p, bounds, Ustate_local, x0, cell);

	while(error >= tolerance) {
		pressureeq(p, Ustate_local, pfunc, x0, cell);
		delta = pfunc[0] / pfunc[1];
		if(fabs(delta) <= bounds[1] - bounds[0]) {
			error = fabs(delta);
		} else {
			error = bounds[1] - bounds[0];
		}

		if(p - delta < bounds[0] || p - delta > bounds[1] || isnan(delta)) {
			pressureeq((bounds[1] + bounds[0]) / 2., Ustate_local, pfunc, x0, cell);
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
		if(iterations > 2000) {
			printf("psolve failed %f %i \n", x0, cell[0]);
			break;
		}
	}

	return p;
}


void cons_to_prim(double *Ustate_ptr, double *Pstate_ptr, double x0) {
	int cell[3];
	double *Ustate_local;
	double *Pstate_local;
	double D, tau, p, W, v_sum, u0;
	double v_i[3];
	double v_i_upper[3];
	double *shift_vect;
	int p_address;	
	for(int i = ghost_num; i < x1cellnum + ghost_num; i++) {
		cell[0] = i;
		for(int j = ghost_num; j < x2cellnum + ghost_num; j++) {
			cell[1] = j;
			for(int k = ghost_num; k < x3cellnum + ghost_num; k++) {
				cell[2] = k;

				v_sum = 0.;
			
				Ustate_local = Ustate_ptr + U_offset(cell, 0);
				Pstate_local = Pstate_ptr + P_offset(cell, 0);

				p = psolve(Ustate_local, Pstate_local, x0, cell);
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

				W = pow(1.0 - v_sum, -.5);
				u0 = W / lapse(x0, cell);
				shift_vect = shift(x0, cell);

				p_address = P_offset(cell, 0);
				*(Pstate_ptr + p_address) = D / W;
				*(Pstate_ptr + p_address + 1) = p;
				*(Pstate_ptr + p_address + 2) = u0;
				*(Pstate_ptr + p_address + 3) = u0 * (lapse(x0, cell) * v_i_upper[0] - shift_vect[0]);
				*(Pstate_ptr + p_address + 4) = u0 * (lapse(x0, cell) * v_i_upper[1] - shift_vect[1]);
				*(Pstate_ptr + p_address + 5) = u0 * (lapse(x0, cell) * v_i_upper[2] - shift_vect[2]);
				
				free(shift_vect);
			}		
		}		
	}
}

