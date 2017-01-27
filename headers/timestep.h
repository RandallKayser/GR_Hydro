void timestep_basic(double *Ustate_ptr, double *Pstate_ptr, double dt, double x0);

void RK2(double *Ustate_ptr, double *Pstate_ptr, double dt);
void RK3(double *Ustate_ptr, double *Pstate_ptr, double dt);
void RK3_TVD(double *Ustate_ptr, double *Pstate_ptr, double dt);
void RK4(double *Ustate_ptr, double *Pstate_ptr, double dt);