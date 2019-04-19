#ifndef __StiefelScheiffele_h__
#define __StiefelScheiffele_h__

#include <dace/dace.h>
using namespace std;
using namespace DACE;

unsigned int funcs = 0;

// initial conditions
namespace SS {
    // km, km/s
    AlgebraicVector<double> x0 = { 0.000000, -5888.972700, -3400.00000, 10.691338, 0.000000, 0.000000 };
    const double epsR = 1e-3, epsV = 1e-4;
    double t0 = 0.;
    //double tf = 288.12768941 * 24. * 3600.;
    //double tf = 288.12768941 * 24. * 3600. / 100;   // Shortened version
    double tf = 288.12768941 * 24. * 3600. / 50;   // Less shortened version
    //AlgebraicVector<double> Ref_Sol ({-24219.05011593605201960070788, 227962.10637302200887202306088, 129753.44240008247047344318589});
    //AlgebraicVector<double> Ref_Sol ({-0.006472054691587e5, 2.290183472899976e5, 1.322860577190283e5});     // Shortened version
    AlgebraicVector<double> Ref_Sol ({-672.68534351214767, -6099.8133178135287, -3522.5764770819265});     // Less shortened version (only computed with tol=1e-13)
}

template <class T>
AlgebraicVector<T> StiefelScheiffele(AlgebraicVector<T> x, double t)
{
    funcs++;
    AlgebraicVector<T> deriv(6);
    AlgebraicVector<T> ap(3);

    AlgebraicVector<T> R = x.extract(0,2);
    AlgebraicVector<T> V = x.extract(3,5);

    // J2 Perturbation
    T Rho = vnorm(R);

    const double R_p = 6371.2200;
    const double J2 = 0.00108265;
    const double GM = 3.98601e5;

    T cost_J2 = -1.5 * GM * J2 * R_p*R_p / pow(Rho, 7);
    T Rho2 = Rho * Rho;
    T z2 = R[2] * R[2];

    ap[0] = cost_J2 * R[0] * (    Rho2 - 5 * z2);
    ap[1] = cost_J2 * R[1] * (    Rho2 - 5 * z2);
    ap[2] = cost_J2 * R[2] * (3 * Rho2 - 5 * z2);

    // Lunar Perturbation
    const double mu_m = 4.90266e3;
    const double d_m = 384400;
    const double w_m = 2.665315780887e-6;

    AlgebraicVector<T> RM(3);
    AlgebraicVector<T> RMrel(3);

    RM[0] =  d_m * sin(w_m * t);
    RM[1] = -d_m * 0.5 * sqrt(3.) * cos(w_m * t);
    RM[2] = -d_m * 0.5 * cos(w_m * t);

    RMrel = R - RM;

    T Rho_M = vnorm(RMrel);
    T Rho_M3 = Rho_M * Rho_M * Rho_M;
    const double d_m3 = d_m * d_m * d_m;

    ap[0] += -mu_m * (RMrel[0] / Rho_M3 + RM[0] / d_m3);
    ap[1] += -mu_m * (RMrel[1] / Rho_M3 + RM[1] / d_m3);
    ap[2] += -mu_m * (RMrel[2] / Rho_M3 + RM[2] / d_m3);

    const double mu = GM;

    T aux = mu / (Rho*Rho*Rho);

    deriv[0] =  x[3];
    deriv[1] =  x[4];
    deriv[2] =  x[5];
    deriv[3] = -aux * x[0] + ap[0];
    deriv[4] = -aux * x[1] + ap[1];
    deriv[5] = -aux * x[2] + ap[2];

    return deriv;
}
#endif