#ifndef INCLUDE_CONSTANTS
#define INCLUDE_CONSTANTS 1

const static int x1cellnum = 250;
const static int x2cellnum = 1;
const static int x3cellnum = 1;
const static int x1type = 2;
const static int x2type = 0;
const static int x3type = 0;

const static int DIM_NUM = 4;

const static double t_init = 0.;
const static double t_final = 100.;

const static double gamma_const = 5./3.;
const static double flux_limiter_theta = 1.5;
const static double cfl_number = .2;
const static int ghost_num = 1;

#endif