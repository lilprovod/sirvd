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

void euler_step(
    double      t,
    double      h,
    double*     x,
    int         dim,
    RhsFunc     f,
    const void* params
)
{
    double dxdt[dim];

    f(t, x, dxdt, params);

    for (int i = 0; i < dim; i++) {
        x[i] += h * dxdt[i];
    }
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
    /*
        k1 - скорость системы в начальный момент t
        k2 - скорость системы через шаг h/2 (в момент t + h/2)

        x_middle - состояние системы посередине (t + h/2); рассчитывается через начальную скорость k1

        x (xn + 1, итоговый) - состояние системы через шаг (t + h); рассчитывается через серединную скорость k2 (прогноз)
    */

    double k1[dim];
    double k2[dim];

    double x_middle[dim];

    f(t, x, k1, params);        // k1 = f(tn, xn)

    for (int i = 0; i < dim; i++) {
        x_middle[i] = x[i] + (double)(h / 2.0) * k1[i];     // x_mid = xn + h/2 * k1
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
    /*
        k1 - скорость в начале (t)          k1 = f(t,       x)
        k2 - скорость в середине (t + h/2)  k2 = f(t + h/2, x + (h/2)*k1)
        k3 - скорость в конце (t + h)       k3 = f(t + h,   x - h*k1 + 2*h*k2)

        x (xn + 1, итоговый) - состояние системы через t + h (прогноз)
        формула: xn + 1 = xn + h/6 (k1 + 4k2 + k3)
    */

    double k1[dim];
    double k2[dim];
    double k3[dim];

    double x_middle[dim];
    double x_end[dim];

    f(t, x, k1, params);

    for (int i = 0; i < dim; i++) {
        x_middle[i] = x[i] + (double)(h / 2.0) * k1[i];
    }

    f(t + h / 2.0, x_middle, k2, params);

    for (int i = 0; i < dim; i++) {
        x_end[i] = x[i] - h * k1[i] + (2.0 * h) * k2[i];
    }

    f(t + h, x_end, k3, params);

    for (int i = 0; i < dim; i++) {
        x[i] = x[i] + (double)(h / 6.0) * (k1[i] + 4.0 * k2[i] + k3[i]);
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
    /*
        k1 - скорость в начале (t)                      k1 = f(t,       x)
        k2 - скорость в середине через k1 (t + h/2)     k2 = f(t + h/2, x + (h/2) * k1)
        k3 - скорость в середине через k2 (t + h/2)     k3 = f(t + h/2, x + (h/2) * k2)
        k4 - скорость в конце (t + h)                   k4 = f(t + h,   x + h * k3)

        x (xn + 1, итоговый) -- состояние системы через t + h (прогнозируемое)
        xn + 1 = xn + h/6 * (k1 + 2k2 + 2k3 + k4)
    */

    double k1[dim];
    double k2[dim];
    double k3[dim];
    double k4[dim];

    double x_mid_k1[dim];
    double x_mid_k2[dim];
    double x_end[dim];

    f(t, x, k1, params);
    
    for (int i = 0; i < dim; i++) {
        x_mid_k1[i] = x[i] + (double)(h / 2.0) * k1[i];
    }

    f(t + h/2.0, x_mid_k1, k2, params);

    for (int i = 0; i < dim; i++) {
        x_mid_k2[i] = x[i] + (double)(h / 2.0) * k2[i];
    }

    f(t + h/2.0, x_mid_k2, k3, params);

    for (int i = 0; i < dim; i++) {
        x_end[i] = x[i] + h * k3[i];
    }

    f(t + h, x_end, k4, params);

    for (int i = 0; i < dim; i++) {
        x[i] = x[i] + (double)(h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
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
    double*     x6
)
{
    /**
     * @param[out] x5
     * @param[out] x6
     */

    enum { S = 9 }; // кол-во стадий

    static const double c[S] = {
        0.0,
        9.0  / 50.0,
        1.0  / 6.0,
        1.0  / 4.0,
        53.0 / 100.0,
        3.0  / 5.0,
        4.0  / 5.0,
        1.0,
        1.0
    };

    static const double a[S][S] = {
        { 0 },
        { 9.0 / 50.0, 0 },
        { 29.0/324.0, 25.0/324.0, 0 },
        { 1.0/16.0, 0.0, 3.0/16.0, 0 },
        { 79129.0/250000.0, 0.0, -261237.0/250000.0, 19663.0/15625.0, 0 },
        { 1336883.0/4909125.0, 0.0, -25476.0/30875.0, 194159.0/185250.0,
          8225.0/78546.0, 0 },
        { -2459386.0/14727375.0, 0.0, 19504.0/30875.0, 2377474.0/13615875.0,
          -6157250.0/5773131.0, 902.0/735.0, 0 },
        { 2699.0/7410.0, 0.0, -252.0/1235.0, -1393253.0/3993990.0,
          236875.0/72618.0, -135.0/49.0, 15.0/22.0, 0 },
        { 11.0/144.0, 0.0, 0.0, 256.0/693.0, 0.0, 125.0/504.0, 125.0/528.0,
          5.0/72.0, 0.0 }
    };

    static const double b5[S] = {
        28.0/477.0, 0.0, 0.0, 212.0/441.0, -312500.0/366177.0,
        2125.0/1764.0, 0.0, -2105.0/35532.0, 2995.0/17766.0
    };
    static const double b6[S] = {
        11.0/144.0, 0.0, 0.0, 256.0/693.0, 0.0, 125.0/504.0, 125.0/528.0, 5.0/72.0, 0.0
    };


    double k[S][dim];
    double x_temp[dim];

    // обход по стадиям
    f(t, x, k[0], params);

    for (int i = 1; i < S; i++) {
        // x_t = x + h * sum (a_ij * kj)
        for (int d = 0; d < dim; d++) {
            double acc = 0.0;
            for (int j = 0; j < i; j++) {
                acc += a[i][j] * k[j][d];
            }

            x_temp[d] = x[d] + h * acc;
        }

        f(t + c[i] * h, x_temp, k[i], params);
    }

    for (int d = 0; d < dim; d++) {
        double sum5 = 0.0;
        double sum6 = 0.0;

        for (int i = 0; i < S; i++) {
            sum5 += b5[i] * k[i][d];
            sum6 += b6[i] * k[i][d];
        }

        x5[d] = x[d] + h * sum5;
        x6[d] = x[d] + h * sum6;
    }
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
    /*
        строим два решения порядка 5 и 6, используем промежуточные k

        оценка ошибки e = x_6 - x_5
        err = max(|e_i| / (atol + rtol * max(|x_i|, |x_6_i|))) для всех i

        если err <= tol ->  x = x_6, t = t + h;
        если err > tol ->   уменьшаем шаг h = h * safety * (tol / err)^(1/p+1) и повторить шаг
    */

    const double safety     = 0.9;
    const double min_factor = 0.2;
    const double max_factor = 5.0;

    // контроль шага
    double hh = *h;
    if (*t + hh > t_end) {
        hh = t_end - *t;
    }
    hh = clamp(hh, h_min, h_max);

    double x5[dim];
    double x6[dim];

    // X5, X6 решения
    rk56_try_step(*t, hh, x, dim, f, params, x5, x6);

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
        return 1;
    }

    *h = clamp(hh * factor, h_min, h_max);
    return 0;
}