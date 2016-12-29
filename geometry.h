#include "constants.h"

double *get_position(int cell[DIM_NUM-1]);
double metric_full(int dir1, int dir2, double x0, int cell[DIM_NUM-1]);
double inner_product(double a[DIM_NUM], double b[DIM_NUM], double x0, int cell[DIM_NUM-1]);
double metric_det(int dir1, int dir2, double x0, int cell[DIM_NUM-1]);

double metric_space(int dir1, int dir2, double x0, int cell[DIM_NUM-1]);
double lapse(double x0, int cell[DIM_NUM-1]);
double *shift(double x0, int cell[DIM_NUM-1]);
double inner_product_space(double a[DIM_NUM-1], double b[DIM_NUM-1]);

double lower_index_space(double v[DIM_NUM-1], int index);

double volume(double x0, int cell[DIM_NUM-1);
