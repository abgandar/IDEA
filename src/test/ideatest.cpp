#include "idea.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace DACE;

template <class T>
class myRHS : public RHS<T>
{
  public:
    T lambda;
    myRHS(T l) : lambda(l) { };
    DACE::AlgebraicVector<T> f(double t, const DACE::AlgebraicVector<T> &x);
};

template <class T>
DACE::AlgebraicVector<T> myRHS<T>::f(double t, const DACE::AlgebraicVector<T> &x)
{
    return DACE::AlgebraicVector<T>({x[1], -x[0]}) * vnorm(x) * lambda;
}

template <class T>
class Kepler : public RHS2<T>
{
  public:
    T mu;
    Kepler(T m) : mu(m){};
    DACE::AlgebraicVector<T> f(double t, const DACE::AlgebraicVector<T> &x);
};

template <class T>
DACE::AlgebraicVector<T> Kepler<T>::f(double t, const DACE::AlgebraicVector<T> &x)
{
    DACE::AlgebraicVector<T> res = x.extract(0,x.size()/2-1);
    T r = vnorm(res);
    res = -res * mu / pow(r,3);
    return res;
}

void Draw2D(AlgebraicVector<DA> v, int N = 100)
{
    AlgebraicVector<double> a(N), bp(N, 1.), bm(N, -1.), x, y;
    for (int i = 0; i < N; i++)
    {
        a[i] = -1. + i * 2. / (N - 1);
    }

    x.insert(x.end(), a.begin(), a.end());
    y.insert(y.end(), bp.begin(), bp.end());

    x.insert(x.end(), bp.begin(), bp.end());
    a = -a;
    y.insert(y.end(), a.begin(), a.end());

    x.insert(x.end(), a.begin(), a.end());
    y.insert(y.end(), bm.begin(), bm.end());
    a = -a;

    x.insert(x.end(), bm.begin(), bm.end());
    y.insert(y.end(), a.begin(), a.end());

    auto r = v.eval<AlgebraicVector<double>>({x, y});

    for (int n = 0; n < r[0].size(); n++)
        cout << r[0][n] << "   " << r[1][n] << endl;
    cout << endl;
}


int main(void)
{
    DA::init(5, 2); // initialize DACE for 20th-order computations in 1 variable

    myRHS<DA> rhsDA(1.);
    myRHS<double> rhsDouble(1.);

    Kepler<DA> KeplerDA(1.);
    Kepler<double> KeplerDouble(1.);

    double eps = 1e-2;
    double h = 0.001;
    double t0 = 0;
    double tf = M_PI * 5;

    AlgebraicVector<DA> x0({1 + eps * DA(1), 1 + eps * DA(2)}); // initialize x as DA
    AlgebraicVector<DA> x;

    AlgebraicVector<double> y0 = x0.eval<double>({1, 0.5});
    AlgebraicVector<double> y;

    //x = Euler(t0, tf, h, x0, rhsDA);
    //y = Euler(t0, tf, h, y0, rhsDouble);


    AlgebraicVector<DA> rvDA({1 + eps*DA(1), eps*DA(2), 0, 1.2});
    //AlgebraicVector<double> rvdouble = rvDA.eval<double> ({1, 0.5});

    x = rvDA;
    int n = 15;
    double t = t0;
    for (int i=0; i<n; i++){
        x = Euler(t, t+(tf-t0)/n, h, x, KeplerDA);
        t = t + (tf-t0)/n;
        Draw2D(x);
    }
    //y = Euler(t0, tf, h, rvdouble, KeplerDouble);

    //cout << tf << endl
    //     << x0 << endl
    //     << x << endl;

    //cout << vnorm(y - x.eval<double>({1, 0.5})) << endl;

}