#include <dace/dace.h>
using namespace std;
using namespace DACE;


// perform a single predict/correct loop and output the new B values
// y, F, G are temporary storage passed to avoid reallocation. y[0], F[0] must be state and RHS at initial time
// c is the array of precomputed coefficients
template <class T, size_t NB, size_t NH>
AlgebraicVector<T> ER_predict_correct( AlgebraicVector<T> &B, double t, double dt, AlgebraicVector<T> &y, AlgebraicVector<T> &F, AlgebraicVector<T> &G, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), double (&c)[NB][NB], double (&h)[NH] )
{
    // predictor
    for( unsigned int i = 1; i < NH; i++ )
        y[i] = y[0] + h[i]*dt*(F[0] + h[i]*(B[0]/2.0 + h[i]*(B[1]/3.0 + h[i]*(B[2]/4.0 + h[i]*B[3]/5.0))));
    // corrector
    for( unsigned int i = 1; i < NH; i++ )
        F[i] = RHS(y[i], t+h[i]*dt);
    G[0] = (F[1]-F[0])/(h[1]-h[0]);
    G[1] = ((F[2]-F[0])/(h[2]-h[0])-G[0])/(h[2]-h[1]);
    G[2] = (((F[3]-F[0])/(h[3]-h[0])-G[0])/(h[3]-h[1])-G[1])/(h[3]-h[2]);
    G[3] = ((((F[4]-F[0])/(h[4]-h[0])-G[0])/(h[4]-h[1])-G[1])/(h[4]-h[2])-G[2])/(h[4]-h[3]);
    AlgebraicVector<T> BB(NB);
    BB[0] = c[0][0]*G[0] + c[1][0]*G[1] + c[2][0]*G[2] + c[3][0]*G[3];
    BB[1] =                c[1][1]*G[1] + c[2][1]*G[2] + c[3][1]*G[3];
    BB[2] =                               c[2][2]*G[2] + c[3][2]*G[3];
    BB[3] =                                              c[3][3]*G[3];
    return BB;
}

template <class T>
AlgebraicVector<T> DAEverhartRadau(double t0, double t1, double hh, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*DARHS)(AlgebraicVector<T>, double), AlgebraicVector<double> (*RHS)(AlgebraicVector<double>, double), vector<AlgebraicVector<T>> &res, bool save)
{
    // Gauss Radau spacings
    const double h[5] = {0.0, 0.13975986434378055, 0.41640956763108318, 0.72315698636187617, 0.94289580388548232};
    const unsigned int NH = sizeof(h)/sizeof(h[0]);

    double t = t0;
    AlgebraicVector<T> B[4], X = x0;
    AlgebraicVector<double> dB[4], dX = cons(x0);
    const unsigned int NB = sizeof(B)/sizeof(B[0]);

    for( unsigned int i = 0; i < NB; i++ )
    {
        B[i] = 0.0*X;
        dB[i] = 0.0*dX;
    }

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
    AlgebraicVector<T> y[NH];
    AlgebraicVector<T> F[NH];
    AlgebraicVector<T> G[NB];

    if( save )
        res.push_back( X );

    while( t < t1*(1.0-1e-14) )
    {
        steps++;
        const double dt = min(t1-t, hh);

        y[0] = X;
        F[0] = RHS(X, t);
        double db = 1.0;
        unsigned int n = 0;
        for( n = 0; n < (steps < 3 ? 20 : 20) && db > 1e-14; n++ )
        {
            const AlgebraicVector<T> B_old = B[3];
            B = ER_predict_correct(B, t, dt, y, F, G, RHS, c, h);
            db = vnorm(B_old[3]-B[3])/vnorm(F[0]);
        }
        X += dt*(F[0] + B[0]/2.0 + B[1]/3.0 + B[2]/4.0 + B[3]/5.0);
        t += dt;

        // update next B guesses
        const double Q = 1.0, Q2 = Q*Q, Q3 = Q2*Q, Q4 = Q3*Q;
        B[0] = Q *(B[0] + 2.0*B[1] + 3.0*B[2] + 4.0*B[3]);
        B[1] = Q2*(           B[1] + 3.0*B[2] + 6.0*B[3]);
        B[2] = Q3*(                      B[2] + 4.0*B[3]);
        B[3] = Q4*(                                 B[3]);

        // save step
        if( save )
            res.push_back( X );

        const double progress = 100.00 * (t) / (t1-t0);
        if( steps % 1000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return X;
};
