/*RK Methods*/

#ifndef RK_H
#define RK_H

#define H_MIN 1e-6
#define H_MAX 10.0


typedef void (*RhsFunc)(double t, const double* x,
                        double* dxdt, const void* params);

double clamp(double v, double lo, double hi);

void euler_step(double t, double h,
                double* x, int dim,
                RhsFunc f, const void* params);

void rk2_step(double t, double h,
              double* x, int dim,
              RhsFunc f, const void* params);

void rk3_step(double t, double h,
              double* x, int dim,
              RhsFunc f, const void* params);

void rk4_step(double t, double h,
              double* x, int dim,
              RhsFunc f, const void* params);

void rk56_try_step(
    double t, double h,
    double* x, int dim,
    RhsFunc f, const void* params,
    double* x5, double* x6);

int rk56_adaptive_step(
    double* t, double* h,
    double h_min, double h_max,
    double* x, int dim,
    RhsFunc f, const void* params,
    double t_end,
    double atol, double rtol);

#endif