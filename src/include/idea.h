#include <dace/dace.h>

template <class T>
class RHS
{
  public:
    virtual DACE::AlgebraicVector<T> operator()(double t, const DACE::AlgebraicVector<T> &x);
    virtual DACE::AlgebraicVector<T> f(double t, const DACE::AlgebraicVector<T> &x) = 0;
};

template <class T>
DACE::AlgebraicVector<T> RHS<T>::operator()(double t, const DACE::AlgebraicVector<T> &x)
{
    return this->f(t, x);
}

template <class T>
class RHS2: public RHS<T>
{
  public:
    DACE::AlgebraicVector<T> operator()(double t, const DACE::AlgebraicVector<T> &x);
};

template <class T>
DACE::AlgebraicVector<T> RHS2<T>::operator()(double t, const DACE::AlgebraicVector<T> &x)
{
    DACE::AlgebraicVector<T> res(x.size()), rhs = this->f(t, x);
    const unsigned int DIM = x.size()/2;
    for(unsigned int i = 0; i < DIM; i++)
    {
        res[i] = x[i+DIM];
        res[i+DIM] = rhs[i];
    }
    return res;
}

template <class T>
DACE::AlgebraicVector<T> Euler(double t0, double tf, double h, const DACE::AlgebraicVector<T> &X0, RHS<T> &rhs)
{
    DACE::AlgebraicVector<T> X = X0;
    int steps = (tf - t0) / h;

    for (int i = 0; i < steps; i += 1)
    {
        X = X + h * rhs(t0 + h * i, X);
    }

    return X;
};


