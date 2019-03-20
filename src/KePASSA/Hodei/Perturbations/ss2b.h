#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace DACE;
using namespace std;

template <class T>
AlgebraicVector<T> SS2B_perturbation (AlgebraicVector<T> &R, AlgebraicVector<T> &V, T Time)
{
    // J2 Perturbation
    T Rho = vnorm(R);
    
    const double R_p = 6371.2200;
    const double J2 = 0.00108265;
    const double GM = 3.98601e5;

    T cost_J2 = -1.5 * GM * J2 * R_p * R_p / pow(Rho, 7);
    T Rho2 = Rho * Rho;
    T z2 = R[2] * R[2];

    AlgebraicVector<T> Contribution(3);
    AlgebraicVector<T> ap(3);

    Contribution[0] = cost_J2 * R[0] * (    Rho2 - 5 * z2);
    Contribution[1] = cost_J2 * R[1] * (    Rho2 - 5 * z2);
    Contribution[2] = cost_J2 * R[2] * (3 * Rho2 - 5 * z2);

    ap = ap + Contribution;

    // Lunar Perturbation
    const double mu_m = 4.90266e3;
    const double d_m = 384400;
    const double w_m = 2.665315780887e-6;

    AlgebraicVector<T> RM(3);
    AlgebraicVector<T> RMrel(3);

    RM[0] =  d_m * sin(w_m * Time);
    RM[1] = -d_m * 0.5 * sqrt(3.) * cos(w_m * Time);
    RM[2] = -d_m * 0.5 * cos(w_m * Time);

    RMrel = R - RM;
    
    T Rho_M = vnorm(RMrel);
    T Rho_M3 = Rho_M * Rho_M * Rho_M;
    double d_m3 = d_m * d_m * d_m;

    Contribution[0] = -mu_m * (RMrel[0] / Rho_M3 + RM[0] / d_m3);
    Contribution[1] = -mu_m * (RMrel[1] / Rho_M3 + RM[1] / d_m3);
    Contribution[2] = -mu_m * (RMrel[2] / Rho_M3 + RM[2] / d_m3);

    ap = ap + Contribution;

    return ap;
}
