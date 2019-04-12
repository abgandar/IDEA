#include <dace/dace.h>
#include <cmath>
using namespace std;
using namespace DACE;


template <class T>
AlgebraicVector<T> Elliptic_ODE(AlgebraicVector<T> x, double t)
{
    funcs++;

    const double GM = 1.; // ODE in non-dimensional form

    AlgebraicVector<T> deriv(4);
    AlgebraicVector<T> R = x.extract(0,1);
    
    T Rho = vnorm(R);
    T aux = GM / (Rho*Rho*Rho);
    
    deriv[0] =   x[2];
    deriv[1] =   x[3];
    deriv[2] = - aux * x[0];
    deriv[3] = - aux * x[1];

    return deriv;
}



void Elliptic_IC(double rp, double ecc, AlgebraicVector<double>& x0, double& T)
{
    const double GM = 1.; // ODE in non-dimensional form

    double a  = rp / (1-ecc);
    double vp = sqrt(GM * (1+ecc) / rp);
    
    // Period of the orbit
    T = 2 * M_PI * sqrt(a*a*a/GM);
    
    // State Vector
    x0[0] = rp;     // X
    x0[1] = 0.;     // Y
    x0[2] = 0.;     // Vx
    x0[3] = vp;     // Vy
}


