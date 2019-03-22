#ifndef __Midpoint_h__
#define __Midpoint_h__

#include <dace/dace.h>
#include <cmath>
using namespace std;
using namespace DACE;

template <class T>
AlgebraicVector<T> Midpoint(double t0, double t1, double h, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), vector<AlgebraicVector<T>> &res, bool save)
{
    const unsigned int steps = round((t1-t0)/h);    // round number of steps to nearest int
    const double dt = (t1-t0)/steps;
    AlgebraicVector<T> z0 = x0, z1 = x0 + dt*RHS(z0, t0);

    if( save )
    {
        res.push_back( z0 );
        res.push_back( z1 );
    }

    for( unsigned int i = 1; i < steps; i++ )
    {
        auto z2 = z0 + 2*dt*RHS(z1, t0 + i*dt);
        z0 = z1;
        z1 = z2;
        if( save )
            res.push_back( z1 );
    }

    z1 = 0.5*( z0 + z1 + dt*RHS(z1, t1));
    if( save )
        res.push_back( z1 );

    return z1;
};

template <class T>
AlgebraicVector<T> Midpoint_RK(double t0, double t1, double h, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), vector<AlgebraicVector<T>> &res, bool save)
{
    const double C[2]    = { 0.0, 0.5 };
    const double B[2]    = { 0.0, 1.0 };
    const double A[2][2] = {{0.,   0.},
                            {0.5,  0.}};
    const unsigned int NS = sizeof(C)/sizeof(C[0]);

    AlgebraicVector<T> Xk = x0, X, K[NS];
    unsigned int steps = 0;
    double t = t0;

    if( save )
        res.push_back( Xk );

    while( t < t1*(1.0-1e-14) )
    {
        steps++;
        const double dt = min(t1-t, h);

        for( unsigned int i = 0; i < NS; i++ )
        {
            X = Xk;
            for( unsigned int j = 0; j < i; j++ )
            {
                if(A[i][j] != 0.0)
                    X += A[i][j]*K[j];
            }
            K[i] = dt*RHS(X, t + C[i]*dt);
        }

        // Estimate y(x+dt) of the highest order available
        for( unsigned int i = 0; i < NS; i++ )
        {
            Xk += B[i] * K[i];
        }

        // Update value of independent variables
        t += dt;

        // save step
        if( save )
            res.push_back( Xk );

        const double progress = 100.00 * (t) / (t1-t0);
        if( steps % 5000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return Xk;
};
#endif