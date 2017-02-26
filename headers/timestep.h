void forward_euler(double *Ustate_ptr, double *U1state_ptr, double *Pstate_ptr, double *x0);
void RK2_TVD(double *Ustate_ptr, double *U1state_ptr, double *Pstate_ptr, double *x0);
void RK3_TVD(double *Ustate_ptr, double *U1state_ptr, double *Pstate_ptr, double *x0);

void set_RK(int type);