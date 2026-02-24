/*SIRVD Model*/

#ifndef SIRVD_H
#define SIRVD_H

#define SIRVD_DIM 5     // размерность системы

typedef struct {
    double N;
    double alpha, beta, gamma, delta, sigma;
} SIRVDParams;

/*
    Вычисляет ОДУ вида dx/dt = f(t, x) для модели SIRVD.

    x: [S, I, R, V, D]
*/
void sirvd_rhs(double t, const double x[SIRVD_DIM],
               double dxdt[SIRVD_DIM],
               const SIRVDParams* p);

void sirvd_rhs_wrap(double t, const double* x,
                    double* dxdt, const void* params);

#endif