#ifndef INCLUDE_CONSTANTS
#define INCLUDE_CONSTANTS 1

const static int x0cellnum = 512;
const static int x1cellnum = 16;
const static int x2cellnum = 16;
const static int x3cellnum = 16;
const static int DIM_NUM = 4;
const static double gamma = 4./3.;
const static double flux_limiter_theta = 1.5;
const static double cfl_number = .2;

const static double mass = 1;
const static double eps = 1e-6;
const static double r_min = 2. * mass + eps;
const static double r_max;

#endif