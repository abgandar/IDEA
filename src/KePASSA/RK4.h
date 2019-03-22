#ifndef __RK4_h__
#define __RK4_h__

#include <dace/dace.h>
using namespace std;
using namespace DACE;

template <class T>
AlgebraicVector<T> RK4(double t0, double t1, double h, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), vector<AlgebraicVector<T>> &res, bool save)
{
    const double C[4]    = { 0.0, 0.5, 0.5, 1.0 };
    const double B[4]    = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
    const double A[4][4] = {{0.,   0.,  0.,  0.0},
                            {0.5,  0.,  0.,  0.0},
                            {0.0,  0.5, 0.,  0.0},
                            {0.0,  0.0, 1.0, 0.0}};
    const unsigned int NS = sizeof(C)/sizeof(C[0]);

    AlgebraicVector<T> Xk = x0, X, K[NS];
    unsigned int steps = 0;
    double t = t0;

    if( save )
        res.push_back( Xk );

    while( t < t1 )
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