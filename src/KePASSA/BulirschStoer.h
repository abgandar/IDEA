#ifndef __BulirschStoer_h__
#define __BulirschStoer_h__

#include <dace/dace.h>
#include "Midpoint.h"
#include "RK4.h"
using namespace std;
using namespace DACE;

template <class T, unsigned int ORDER = 4>
AlgebraicVector<T> BulirschStoer(double t0, double t1, double h, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), vector<AlgebraicVector<T>> &res, bool save)
{
    AlgebraicVector<T> eta[ORDER][ORDER], X = x0;
    double n[ORDER+1] = { 2.0, 4.0, 6.0 };
    unsigned int steps = 0;
    double t = t0;

    // step size multipliers
   for( unsigned int i = 3; i < ORDER; i++ )
       n[i] = 2.0*n[i-2];
    // for( unsigned int i = 0; i < ORDER; i++ )
    //     n[i] = 2*i+2;

    if( save )
        res.push_back( X );

    while( t < t1 )
    {
        steps++;
        const double dt = min(t1-t, h);

        // Compute intermediate results
        for( unsigned int i = 0; i < ORDER; i++ )
            eta[i][0] = Midpoint(t, t+dt, dt/n[i], tol, X, RHS, res, false);

        // Build interpolation tables
        for( unsigned int k = 1; k < ORDER; k++ )
            for( unsigned int j = k; j < ORDER; j++ )
                eta[j][k] = eta[j][k-1] + (eta[j][k-1] - eta[j-1][k-1])/(n[j]*n[j]/(n[j-k]*n[j-k]) - 1.0);

        // Estimate y(x+dt) of the highest order available
        X = eta[ORDER-1][ORDER-1];

        // Update value of independent variables
        t += dt;

        // save step
        if( save )
            res.push_back( X );

        const double progress = 100.00 * (t) / (t1-t0);
        if( steps % 5000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return X;
};
#endif