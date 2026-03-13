/*RK Methods*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rk.h"

double clamp(double v, double lo, double hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;

    return v;
}

void rk2_step(
    double      t,
    double      h,
    double*     x,
    int         dim,
    RhsFunc     f,
    const void* params
)
{
    double k1[dim], k2[dim];
    double x_middle[dim];

    f(t, x, k1, params);        // k1 = f(tn, xn)
    for (int i = 0; i < dim; i++) {
        x_middle[i] = x[i] + (h / 2.0) * k1[i];     // x_mid = xn + h/2 * k1
    }

    f(t + h / 2.0, x_middle, k2, params);   // k2 = f(t + h/2, x)
    for (int i = 0; i < dim; i++) {
        x[i] = x[i] + h * k2[i];        // xn+1 = xn + h*k2
    }
}


void rk3_step(
    double      t,
    double      h,
    double*     x,
    int         dim,
    RhsFunc     f,
    const void* params
)
{
    double k1[dim], k2[dim], k3[dim];
    double x_middle[dim], x_end[dim];

    f(t, x, k1, params);
    for (int i = 0; i < dim; i++) {
        x_middle[i] = x[i] + (h / 2.0) * k1[i];
    }

    f(t + h / 2.0, x_middle, k2, params);
    for (int i = 0; i < dim; i++) {
        x_end[i] = x[i] - h * k1[i] + (2.0 * h) * k2[i];
    }

    f(t + h, x_end, k3, params);
    for (int i = 0; i < dim; i++) {
        x[i] = x[i] + (h / 6.0) * (k1[i] + 4.0 * k2[i] + k3[i]);
    }
}


void rk4_step(
    double      t,
    double      h,
    double*     x,
    int         dim,
    RhsFunc     f,
    const void* params
)
{
    double k1[dim], k2[dim], k3[dim], k4[dim];
    double x_mid_k1[dim], x_mid_k2[dim], x_end[dim];

    f(t, x, k1, params);
    for (int i = 0; i < dim; i++) {
        x_mid_k1[i] = x[i] + (h / 2.0) * k1[i];
    }

    f(t + h/2.0, x_mid_k1, k2, params);
    for (int i = 0; i < dim; i++) {
        x_mid_k2[i] = x[i] + (h / 2.0) * k2[i];
    }

    f(t + h/2.0, x_mid_k2, k3, params);
    for (int i = 0; i < dim; i++) {
        x_end[i] = x[i] + h * k3[i];
    }

    f(t + h, x_end, k4, params);
    for (int i = 0; i < dim; i++) {
        x[i] = x[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

void rk56_try_step(
    double      t,
    double      h,
    double*     x,
    int         dim,
    RhsFunc     f,
    const void* params,
    double*     x5,
    double*     x6,
    double*     k1,
    int         use_fsal,
    double*     k_last
)
{
    double k[RK56_S][dim];
    double x_temp[dim];

    if (use_fsal) {
        memcpy(k[0], k1, dim * sizeof(double));
    } else {
        f(t, x, k[0], params);
    }

    for (int i = 1; i < RK56_S; i++) {
        for (int d = 0; d < dim; d++) {
            double acc = 0.0;
            for (int j = 0; j < i; j++) {
                acc += RK56_A[i][j] * k[j][d];
            }

            // x_t = x + h * sum (a_ij * kj)
            x_temp[d] = x[d] + h * acc;
        }

        f(t + RK56_C[i] * h, x_temp, k[i], params); // k_i = f(t + c_i * h, x_i)
    }

    for (int d = 0; d < dim; d++) {
        double sum5 = 0.0;
        double sum6 = 0.0;

        for (int i = 0; i < RK56_S; i++) {
            sum5 += RK56_B5[i] * k[i][d];
            sum6 += RK56_B6[i] * k[i][d];
        }

        x5[d] = x[d] + h * sum5;
        x6[d] = x[d] + h * sum6;
    }

    memcpy(k_last, k[RK56_S - 1], dim * sizeof(double));
}

int rk56_adaptive_step(
    double*     t,
    double*     h,
    double      h_min,
    double      h_max,
    double*     x,
    int         dim,
    RhsFunc     f,
    const void* params,
    double      t_end,
    double      atol,
    double      rtol
)
{
    const double safety     = 0.9;
    const double min_factor = 0.2;
    const double max_factor = 2.0;

    static double k_fsal[64];   // SYSTEM MAX_DIM = 64
    static int use_fsal = 0;

    // контроль шага
    double hh = *h;
    if (*t + hh > t_end) {
        hh = t_end - *t;
    }
    hh = clamp(hh, h_min, h_max);

    double x5[dim], x6[dim];

    // X5, X6 решения
    double k_last[dim];
    rk56_try_step(*t, hh, x, dim, f, params, x5, x6, k_fsal, use_fsal, k_last);

    double err = 0.0;
    for (int i = 0; i < dim; i++) {
        double scale = atol + rtol * fmax( fabs(x[i]), fabs(x6[i]) );
        if (scale == 0.0) scale = 1.0;

        double ei = fabs(x6[i] - x5[i]) / scale;

        if (ei > err) err = ei;
    }
    
    double factor;
    if (err == 0.0) {
        factor = max_factor;
    } else {
        factor = safety * pow(1.0 / err, 1.0 / 6.0);
    }

    factor = clamp(factor, min_factor, max_factor);

    if (err <= 1.0) {
        memcpy(x, x6, (size_t)dim * sizeof(double));
        *t += hh;
        *h = clamp(hh * factor, h_min, h_max);

        memcpy(k_fsal, k_last, dim * sizeof(double));
        use_fsal = 1;

        return 1;
    }

    *h = clamp(hh * factor, h_min, h_max);
    use_fsal = 0;

    return 0;
}