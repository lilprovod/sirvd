/*SIRVD Model*/

#include "sirvd.h"

void sirvd_rhs(
    double             t,
    const double       x[SIRVD_DIM],
    double             dxdt[SIRVD_DIM],
    const SIRVDParams* p
)
{
    (void)t;

    // параметры системы
    const double S = x[0];
    const double I = x[1];
    const double R = x[2];
    const double V = x[3];
    const double D = x[4];

    // население и коэф. системы
    const double N = p->N;
    
    const double alpha = p->alpha;      // vacine coef
    const double beta  = p->beta;       // infect coef
    const double gamma = p->gamma;      // recover coef
    const double delta = p->delta;      // death coef
    const double sigma = p->sigma;      // suscept coef

    dxdt[0] = -(beta * I * S) / N  +  sigma * R  - alpha * S;   // dS/dt
    dxdt[1] =  (beta * I * S) / N  - gamma * I  - delta * I;    // dI/dt
    dxdt[2] = gamma * I - sigma * R;                            // dR/dt
    dxdt[3] = alpha * S;                                        // dV/dt
    dxdt[4] = delta * I;                                        // dD/dt

    return;
}


void sirvd_rhs_wrap(
    double          t,
    const double*   x,
    double*         dxdt,
    const void*     params)
{
    const SIRVDParams* p = (const SIRVDParams*)params;
    sirvd_rhs(t, x, dxdt, p);
}