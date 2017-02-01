#ifndef INCLUDE_CONSTANTS
#define INCLUDE_CONSTANTS 1

const static int x1cellnum = 50;
const static int x2cellnum = 50;
const static int x3cellnum = 3;
const static int x1type = 0;
const static int x2type = 0;
const static int x3type = 0;

const static int DIM_NUM = 4;

const static double t_init = 0.;
const static double t_final = 3e4;

const static double gamma_const = 4./3.;
const static double flux_limiter_theta = 1.5;
const static double cfl_number = .2;
const static int ghost_num = 1;
const static double mass = 1470.;
const static double eps = 1e-6;
const static double r_min = 2. * mass + eps;
const static double r_max = 200 * r_min;

#endif