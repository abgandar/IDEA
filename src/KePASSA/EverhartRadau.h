#include <dace/dace.h>
using namespace std;
using namespace DACE;


template <class T>
AlgebraicVector<T> EverhartRadau(double t0, double t1, double hh, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), vector<AlgebraicVector<T>> &res, bool save)
{
    // Gauss Radau spacings
    const double h[5] = {0.0, 0.13975986434378055, 0.41640956763108318, 0.72315698636187617, 0.94289580388548232};
    const unsigned int NH = sizeof(h)/sizeof(h[0]);

    double t = t0;
    AlgebraicVector<T> B[4], X = x0;
    const unsigned int NB = sizeof(B)/sizeof(B[0]);

    for( unsigned int i = 0; i < NB; i++ )
        B[i] = 0.0*x0;

    double c[NB][NB];
    c[0][0] = 1.0;
    for( unsigned int i = 1; i < NB; i++ )
    {
        c[i][i] = 1.0;
        c[i][0] = -h[i]*c[i-1][0];
        for( unsigned int j = 1; j < i; j++ )
            c[i][j] = c[i-1][j-1] - h[i]*c[i-1][j];
    }

    unsigned int steps = 0;

    if( save )
        res.push_back( X );

    bool first = true;
    while( t < t1*(1.0-1e-14) )
    {
        steps++;
        const double dt = min(t1-t, hh);

        AlgebraicVector<T> y[NH];
        AlgebraicVector<T> F[NH];
        AlgebraicVector<T> G[NB];
        y[0] = X;
        F[0] = RHS(X, t);
        for( unsigned int n = 0; n < (first ? 20 : 4); n++ )
        {
            // predictor
            for( unsigned int i = 1; i < NH; i++ )
                y[i] = X + h[i]*dt*(F[0] + h[i]*(B[0]/2.0 + h[i]*(B[1]/3.0 + h[i]*(B[2]/4.0 + h[i]*B[3]/5.0))));
            // corrector
            for( unsigned int i = 1; i < NH; i++ )
                F[i] = RHS(y[i], t+h[i]*dt);
            G[0] = (F[1]-F[0])/(h[1]-h[0]);
            G[1] = ((F[2]-F[0])/(h[2]-h[0])-G[0])/(h[2]-h[1]);
            G[2] = (((F[3]-F[0])/(h[3]-h[0])-G[0])/(h[3]-h[1])-G[1])/(h[3]-h[2]);
            G[3] = ((((F[4]-F[0])/(h[4]-h[0])-G[0])/(h[4]-h[1])-G[1])/(h[4]-h[2])-G[2])/(h[4]-h[3]);
            B[0] = c[0][0]*G[0] + c[1][0]*G[1] + c[2][0]*G[2] + c[3][0]*G[3];
            B[1] = c[1][1]*G[1] + c[2][1]*G[2] + c[3][1]*G[3];
            B[2] = c[2][2]*G[2] + c[3][2]*G[3];
            B[3] = c[3][3]*G[3];
        }
        X += dt*(F[0] + B[0]/2.0 + B[1]/3.0 + B[2]/4.0 + B[3]/5.0);
        t += dt;
        first = false;

        // save step
        if( save )
            res.push_back( X );

        const double progress = 100.00 * (t) / (t1-t0);
        if( steps % 1000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return X;
};
