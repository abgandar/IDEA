#include <dace/dace.h>
using namespace std;
using namespace DACE;

template <class T>
AlgebraicVector<T> Leapfrog(double t0, double t1, double h, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), vector<AlgebraicVector<T>> &res, bool save)
{
    AlgebraicVector<T> X = x0, xi, vi, ai;
    unsigned int steps = 0;
    double t = t0;

    if( save )
        res.push_back( xi );

    xi = X.extract(0, 2);
    vi = X.extract(3, 5);
    ai = RHS( X, t ).extract(3, 5);
    while( t < t1 )
    {
        steps++;
        const double dt = min(t1-t, h);

        auto xii = xi + vi*dt + 0.5*ai*dt*dt;
        X[0] = xii[0]; X[1] = xii[1]; X[2] = xii[2];
        auto aii = RHS( X, t+dt ).extract(3, 5);
        auto vii = vi + 0.5*(ai + aii)*dt;
        X[3] = vii[0]; X[4] = vii[1]; X[5] = vii[2];

        // Update value of independent variables
        t += dt;
        ai = aii;
        xi = xii;
        vi = vii;

        // save step
        if( save )
            res.push_back( X );

        const double progress = 100.00 * (t) / (t1-t0);
        if( steps % 5000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return X;
};
